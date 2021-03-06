'''
Run a single period model to optimise selection of facilities to open
'''
# import libraries
import os
import yaml
import optimize
import import_location
import matplotlib.pyplot as plt
import numpy as np
import math
import time
import pandas as pd
import geopandas as gpd
from shapely import wkt
import osf
import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

# config_filename = 'denver-set_kpcoef'
config_filename = 'neworleans-set_kpcoef'
# config_filename = 'neworleans-kp_linear_exact'
# config_filename = 'neworleans-pmedian'
# config_filename = 'grid-set_kpcoef'
# config_filename = 'denver-kp_linear_exact'
# config_filename = 'grid-kp_linear_exact'
# config_filename = 'wilmington-set_kpcoef'
# config_filename = 'wilmington-kp_linear_exact'
# config_filename = 'wilmington-pmedian'

def main():
    logger.info('running for {}'.format(config_filename))
    # import the config file
    with open('./config/{}.yml'.format(config_filename)) as file:
        config = yaml.safe_load(file)

    # 8,11, 23,import the location data
    location = import_location.main(config)

    # optimize to identify new facilities to open
    open_optimal = optimize_facility_location(config, location)
    print(open_optimal)

    # plot
    if config['plot']:
        if config['location'] == 'grid':
            plot_grid(config, open_optimal)
        else:
            plot_map(config, open_optimal, location)
    return open_optimal

def optimize_facility_location(config, location):
    # unpack the variables
    origins, destinations = location['origins'], location['destinations']
    populations, distances = location['populations'], location['distances']
    open_current = location['existing']

    #below is code added to adjust the distances in the data frame to be e^(alpha*d[o,d]), which cuts down computation time significantly during optimization
    if config['optimize'] == 'pmedian':
        distances = location['distances']
    else:
        #print('debug100',distances[(371299801001000, 197412)])
        #print('debug101',np.exp(location['alpha']*distances[(371299801001000, 197412)]))
        #distances[(371299801001000, 197412)]=np.exp(location['alpha']*distances[(371299801001000, 197412)])
        for key in distances.keys():
            #distances[key] = np.exp(location['alpha']*distances[key])
            d_key=distances[key]  
            #for distances larger than 63, scip can not handle how big the number e^(alpha*d[o,d]) is, there are only 76 values of d[o,d] larger than 63
            if d_key < 63: distances[key] = np.exp(location['alpha']*distances[key])
            else: 
                #print('key, too big',np.exp(location['alpha']*d_key))
                distances[key]=np.exp(location['alpha']*63)
                #print('too big!!!!',key,np.exp(location['alpha']*63))

    #added to fix error of open total     
    if config['optimize'] == 'pmedian':
        open_total = len(open_current) + config['num_to_open']
    elif config['optimize'] == 'kolmpollak':
        open_total = len(open_current) + config['num_to_open']
    elif config['optimize'] == 'piecewise_linear':
        open_total = len(open_current) + config['num_to_open']                                                                             
    elif config['optimize'] == 'kp_linear_exact':
        open_total = len(open_current) + config['num_to_open']

    # select the optimization algorithm
    if config['optimize'] == 'pmedian':
        open_optimal = optimize.pmedian(origins, destinations, populations,
                                          distances, open_total, open_current)
    elif config['optimize'] == 'kolmpollak':
        open_optimal = optimize.kolmpollak(origins, destinations, populations,
                                             distances, open_total, open_current, location['alpha'])
    elif config['optimize'] == 'piecewise_linear':
        open_optimal = optimize.piecewise_linear(origins, destinations, populations,
                                             distances, open_total, open_current, location['alpha'])                                         
    elif config['optimize'] == 'kp_linear_exact':
        open_optimal = optimize.kp_linear_exact(origins, destinations, populations,
                                             distances, open_total, open_current, location['alpha']) 
    elif config['optimize'] == 'set_kpcoef':
        open_optimal = optimize.set_kpcoef(origins, destinations, populations, 
                                             distances, open_current, location['alpha'], kpcoef=349301)        
                                                                           	
    return(open_optimal)



def plot_grid(config, open_optimal):
    grid_size = config['grid_size']
    fig = plt.figure()
    ax = plt.gca()
    ax.set_ylim(0, grid_size)
    ax.set_xlim(0, grid_size)
    ax = plt.gca().set_aspect("equal")
    plt.xticks(np.arange(0, grid_size + 1, step=1), fontsize=15)
    plt.yticks(np.arange(0, grid_size + 1, step=1), fontsize=15)
    plt.grid(True)
    horiz_grid_coord = list(range(0, grid_size))
    count = 1
    for facility in open_optimal:
        vert_coord = abs(math.floor(facility/grid_size))
        horiz_coord = facility - vert_coord*grid_size
        plt.scatter(horiz_coord + 0.5,
                    abs(math.floor(facility/grid_size) + 0.5))
        plt.annotate("{}".format(count), (horiz_coord + 0.5, abs(math.floor(facility/grid_size) + 0.5)),
                     xytext=(4, 4), textcoords="offset points", ha='center', va='bottom')
        count += 1

    figname = 'fig/{location}-{optimize}'.format(location=config['location'],optimize=config['optimize'])
    plt.savefig(figname + '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.png')



def plot_map(config, open_optimal, location):
    '''
    this currently is only written to work for New Orleans
    '''
    gdf = gpd.read_file('data/new_orleans/raw/tl_2010_22071_tabblock10.shp')
    gdf['origin'] = gdf['STATEFP10'] + gdf['COUNTYFP10'] + gdf['TRACTCE10'] + gdf['BLOCKCE10']
    gdf['origin'] = gdf['origin'].astype('int64')
    # match with origins
    populations = pd.DataFrame.from_dict(location['populations'], orient='index',columns=['population'])
    gdf = gdf.merge(populations, how='right', left_on='origin', right_index=True)
    citymap = gdf.plot(color='white', edgecolor='black')
    xlim = citymap.get_xlim()
    ylim = citymap.get_ylim()
    # plot existing supermarkets
    city_name = config['location'].lower().replace(" ", "_")
    file_path = './data/' + city_name + '/'
    dests = pd.read_csv(file_path + city_name + '-destinations.csv')
    dests['geometry'] = dests['geom'].apply(wkt.loads)
    dests = gpd.GeoDataFrame(dests, crs='epsg:4326')
    existing = dests[dests['dest_type']==config['service']]
    existing.plot(ax=citymap, marker='o',color='blue',markersize=10)
    # new locations
    new = dests[(dests['id_dest'].isin(open_optimal)) & (dests['dest_type']==config['candidate'])]
    new.plot(ax=citymap, marker='x', color='red', markersize=20)
    # map limits
    citymap.set_xlim(xlim)
    citymap.set_ylim(ylim)
    # save
    figname = 'fig/{location}-{optimize}'.format(location=config['location'],optimize=config['optimize'])
    plt.savefig(figname + '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.png')


if __name__ == '__main__':
    start = time.time()
    open_optimal=main()
    elapsed = time.time() - start
    with open('./config/{}.yml'.format(config_filename)) as file:
        config = yaml.safe_load(file)
    location = config['location']
    if location=='grid':
        location = 'grid_' + str(config['grid_size'])
    add = pd.DataFrame([[location, config['optimize'], config['num_to_open'], elapsed,open_optimal]], columns=[
                       'location', 'algorithm', 'number to open', 'computational time', 'optimal stores'])
    output_path = 'computation_time.csv'
    add.to_csv(output_path, mode='a', header=not os.path.exists(output_path))
