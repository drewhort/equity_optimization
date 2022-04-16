import yaml
import main
import pandas as pd
import geopandas as gpd

config_filename = 'new-orleans'
# import config file
with open('./src/query/config/{}.yaml'.format(config_filename)) as file:
    config = yaml.load(file)

# connect to the psql database
db = main.init_db(config)

# download the data for the city
sql = 'SELECT * FROM distxdem;'
dist = pd.read_sql(sql, db["con"])
dist.to_csv('/homedirs/projects/equitable_facility_location/data/{}_population.csv'.format(config_filename))

# download distances
sql = "SELECT * FROM distance;"
dist = pd.read_sql(sql, db['con'])
dist.to_csv('/homedirs/projects/equitable_facility_location/data/{}_distances.csv'.format(config_filename))


# download destinations
sql = "SELECT * FROM destinations"
dest = gpd.GeoDataFrame.from_postgis(sql, db['con'], geom_col='geom')
dest
dest['lon'] = dest.geom.centroid.x
dest['lat'] = dest.geom.centroid.y
dest = dest.rename(columns={'id':'id_dest'})
dest['closed'] = 1
dest.loc[dest.dest_type=='supermarket','closed'] = 0
dest.to_csv('/homedirs/projects/equitable_facility_location/data/{}_destinations.csv'.format(config_filename))




