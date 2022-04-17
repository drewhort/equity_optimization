'''
Import the variables for the location
'''

import math
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import inequalipy as ineq
import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


def main(config):
    if config['location'] == 'grid':
        location = create_grid(config['grid_size'])
    else:
        # import the data from the folder
        city_name = config['location'].lower().replace(" ", "_")
        location = import_files(city_name, config)

    if config['optimize'] == 'kolmpollak':
        # calculate kappa
        location['alpha'] = calculate_kappa(location, config)
    elif config['optimize'] == 'piecewise_linear':
        # calculate kappa
        location['alpha'] = calculate_kappa(location, config)    
    elif config['optimize'] == 'kp_linear_exact':
        # calculate kappa
        location['alpha'] = calculate_kappa(location, config) 
    elif config['optimize'] == 'set_kpcoef':
        # calculate kappa
        location['alpha'] = calculate_kappa(location, config) 

    location['populations'] = location['populations']['population'].to_dict()
    location['distances'] = location['distances']['distance'].to_dict()
    #print('location of distandes',location['distances'])
    return location


def create_grid(grid_size):
    logger.info('generating the grid')
    # create the origins and destinations
    origins = [z for z in range(0, (grid_size*grid_size))]
    destinations = origins.copy()

    # init population weighting
    populations = pd.DataFrame(1, index=origins, columns=['population'])
    populations.index.set_names('origin',inplace=True)

    # calculate distances
    distances = []
    for j in origins:
        for i in destinations:
            # Rounds value down so number of squared in vertical direction can be calculated
            vert_distance = abs(math.floor(i/grid_size) - math.floor(j/grid_size))
            horiz_distance = abs((i % grid_size) - j % grid_size)
            distance = math.sqrt((horiz_distance)**2 + (vert_distance)**2)
            distances.append([j, i, distance])
    distances = pd.DataFrame(distances, columns=['origin','destination','distance'])
    distances.set_index(['origin','destination'], inplace=True)

    return {'origins': origins, 'destinations': destinations, 'populations': populations, 'distances': distances, 'existing': []}


def import_files(city_name, config):
    file_path = './data/' + city_name + '/'

    # populations as a df with origins as index
    logger.info('importing the population')
    populations = pd.read_csv(file_path + city_name + '-population.csv')
    # keep only populated locations
    # populations = populations[populations['U7B001'] > 0]
    # populations.rename(columns={'id_orig': 'origin', 'U7B001':'population'}, inplace=True)
    populations = populations[populations['H7X001'] > 0]
    # print(populations)
    populations.rename(columns={'geoid10': 'origin', 'H7X001':'population'}, inplace=True)
    #what does line below do again? 
    populations = populations[['population', 'origin']]
    # print(populations)
    populations.set_index(['origin'], inplace=True)
   
    # destinations as a list
    logger.info('importing the destinations')
    destinations_df = pd.read_csv(file_path + city_name + '-destinations.csv')
    if config['objective'] == 'restore':
        destinations = destinations_df.loc[destinations_df['dest_type']==config['service'],'id_dest'].unique().tolist()
        # identify open stores
        open_current = destinations_df.loc[(destinations_df['closed'] == 0) & (
            destinations_df['dest_type'] == config['service']),'id_dest'].tolist()
    elif config['objective']=='add':
        destinations = destinations_df.loc[destinations_df['dest_type'].isin([config['service'],config['candidate']]), 'id_dest'].unique().tolist()
        open_current = destinations_df.loc[destinations_df['dest_type']
                                            == config['service'], 'id_dest'].unique().tolist()

    # distances as a df with multi-index: ['origin','destination']
    logger.info('importing the distances')
    distances = pd.read_csv(file_path + city_name + '-distances.csv')
    # get origins that are inhabited
    distances = distances.loc[distances['id_orig'].isin(populations.index)]
    distances = distances.loc[distances['id_dest'].isin(destinations)]
    distances.rename(columns={'id_orig': 'origin', 'id_dest': 'destination'}, inplace=True)
    distances = distances[['origin', 'destination', 'distance']]
    distances.set_index(['origin', 'destination'], inplace=True)
    distances.sort_index(level=[0,1], inplace=True)

    # origins as a list
    origins = distances.index.get_level_values(0).unique().tolist()
    
    logger.info('data imported')
    return {'origins': origins, 'destinations': destinations, 'populations': populations, 'distances': distances, 'existing': open_current}


def calculate_kappa(location, config):
    # get open facilities
    open_current = location['existing']
    # ----> this doesn't work for recovery. in that case, use the ones that were open before disaster.
    if len(open_current) == 0:
        # randomly choose some to open for purpose of calculating kappa
        open_current = random.sample(location['destinations'], max(1, config['grid_size']//2))
    
    # determine the nearest distances
    distances = location['distances']
    #print('distances',distances)
    distances_open = distances.loc[pd.IndexSlice[:, open_current], :]
    nearest = distances_open.groupby('origin')['distance'].min().values
    #print('nearest',nearest)

    # calculate the alpha value
    kappa = ineq.kolmpollak.calc_kappa(nearest, epsilon = config['epsilon'], weights=location['populations'].values)
    alpha = -kappa
    print(alpha)

    return alpha
