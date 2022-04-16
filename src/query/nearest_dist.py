'''
Populate the database for the nearest proximity throughout time
'''
import pandas as pd
import numpy as np
from tqdm import tqdm
import yaml
import main

def determine_nearest(config):
    '''
    determine closest services for time = 0 (initial case), all ids are open
    '''
    db = main.init_db(config)

    con = db['con']
    cursor = con.cursor()

    service = 'supermarket'#config['services']

    sql = """
        SELECT DISTINCT ON (id_orig)
            id_orig, id_dest, distance, duration, dest_type
        FROM
            distance
        WHERE
            id_dest IN (SELECT id FROM destinations WHERE dest_type = %s)
        ORDER BY
            id_orig, distance ASC;
        """
    df = pd.read_sql(sql, db['con'], params = (service,))
    #
    # #get destination ids
    # outs = {}
    #
    # for service in services:
    #     # outs[service] = import_outages(service)
    #     sql = "SELECT id FROM destinations WHERE dest_type = %s;"
    #     dests = pd.read_sql(sql, db['con'], params = (service,))
    #     dest_ids = dests.id.values
    #     outs[service] = {'0':dest_ids}
    #
    # # init the dataframe
    # df = pd.DataFrame(columns = ['id_orig','distance','service', 'time_stamp', 'sim_num', 'metric'])
    # # get the times
    # times = sorted(outs[services[0]].keys()) #times is just ['0'] in this initial case
    # time_stamp = times[0] #because this is the initial case where everything is open
    # # get the distance matrix
    # distances = pd.read_sql('SELECT * FROM distance', db['con'])
    # #making the id_dest column the index
    # distances = distances.set_index('id_dest')
    # #converts distance column from string to float
    # distances.distance = pd.to_numeric(distances.distance)
    #
    # # block ids
    # id_orig = np.unique(distances.id_orig)
    #
    # # loop services
    # for i in tqdm(range(len(services))):
    #     service = services[i]
    #     ids_open = outs[service][time_stamp]
    #     # subset the distance matrix on dest_id
    #     dists_sub = distances.loc[ids_open]
    #     # get the minimum distance
    #     df_min = dists_sub.groupby('id_orig')['distance'].min()
    #     # prepare df to append: This makes distance the name of the column not the series
    #     df_min = df_min.to_frame('distance')
    #     df_min.reset_index(inplace=True)
    #     import code
    #     code.interact(local=locals())
    #
    #     # prepare df to append. Adding these columns as they are needed in simulation.py
    #     df_min['service'] = service
    #     df_min['time_stamp'] = time_stamp
    #     df_min['sim_num'] = 0
    #     df_min['metric'] = 0
    #     # append
    #     df = df.append(df_min, ignore_index=True)
    #sorts by id_orig
    # df.sort_values(by=['id_orig', 'service'], inplace=True)
    # add df to sql, if it exists it will be replaced
    df.to_sql('nearest_dist', con=db['engine'], if_exists='replace', index=False)
    # add index
    # cursor.execute('CREATE INDEX on nearest_dist (time_stamp);')
    # commit
    con.commit()
    con.close()


if __name__ == "__main__":
    config_filename = 'new-orleans'
    # import config file
    with open('./src/query/config/{}.yaml'.format(config_filename)) as file:
        config = yaml.load(file)
    determine_nearest(config)
