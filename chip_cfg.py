import configparser
import pymysql
import sqlite3


def DEPRECATED_read_configuration(path):
    """
        Lee la configuracion del fichero .ini pasado por parametro

        TODO:
            * Hacerlo robusto (frente a errores)
    """
    config = configparser.ConfigParser(interpolation=None)
    config.read(path)


    return config


def root_dir(config, sample, workflow):
    """ Get `rootdir` directory """
    return f"{config.get('DIR', 'analysis_dir')}/{sample}/{workflow}"


# DB

def connect_db(config):
    try:
        db = pymysql.connect(host=config.get('MYSQL', 'host'),
                             user=config.get('MYSQL', 'user'),
                             passwd=config.get('MYSQL', 'passwd'),
                             database=config.get('MYSQL', 'db'),
                             port=config.getint('MYSQL', 'port'))
        return db

    except Exception as e:
        print(f"ERROR: Can't connect to DB : {str(e)}")
        exit(1)


def exec_query(db, query):
    """ Execute qury on MySQL"""
    cursor = db.cursor()
    cursor.execute(query)
    return cursor.fetchall()

def exec_query_to_dict(db, query):
    """ Execute qury on MySQL"""
    cursor = db.cursor(pymysql.cursors.DictCursor)
    cursor.execute(query)
    return cursor.fetchall()


#OLD
def DEPRECATED_root_dir_OLD2(config, workflow, sample, runid):
    """ Get rootdir directory.

        DEPRECATED
    """
    if workflow == 'consensus_reads':
        return f"{config.get('DIR','analysis_dir')}/{sample}/{workflow}/{runid}"
    else:
        exit(1)


def DEPRECATED_exec_query_sqlite(config, query):
    """ Execute query on SQLite"""

    db = config.get('DB', 'db')
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute(query)
    return cur.fetchall()

