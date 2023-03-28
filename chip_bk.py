import shutil
import sys

from TargetPanel import TargetPanel
import chip_cfg
import csv
from io import BytesIO, StringIO
import json
from pathlib import Path
import requests
from cromwell_tools.cromwell_auth import CromwellAuth
from cromwell_tools.cromwell_api import CromwellAPI
import subprocess
import contextlib
import sqlite3
import pymysql
from functools import reduce
import re
# import report as rp
import chip_stats as stats
import pandas as pd

from TargetPanelMySQLConnector import TargetPanelMySQLConnector


# Aux workflow
def is_valid_state(state):
    """ Checks if its a valid state"""
    return state in ('_ready', '_in_progress', '_successful', '_failed', '_aborted')


def get_workflow(step):
    """ Get thw workflow for a given step"""
    switcher = {
        1: 'consensus_reads',
        2: 'variant_discovery',
        3: 'somatic_calling'
    }

    return switcher.get(step)


def get_result_files_suffixes(step):
    """ Returns a list of the result files suffixes for a given step"""

    if step == 1:
        """
        return [".bam", ".bai", ".base_coverage", ".flagstat", \
                     ".hs_metrics", ".markdup.grouped.base_coverage", \
                     ".markdup.grouped.hist", ".markdup.grouped.hs_metrics", \
                     ".markdup.grouped.RX.barcode", \
                     ".markdup.grouped.sorted.consensus.mapped.bam", \
                     ".markdup.grouped.sorted.consensus.mapped.bai", \
                     ".markdup.grouped.target_coverage", \
                     ".markdup.markduplicates_metrics", \
                     ".markdup.umi_metrics", f".target_coverage"]
        """
        return [".bam", ".bai", ".base_coverage", ".target_coverage", ".hs_metrics", ".flagstat",
                ".markdup.grouped.RX.barcode", ".markdup.markduplicates_metrics", ".markdup.umi_metrics",
                ".markdup.grouped.hist", ".markdup.grouped.base_coverage", ".markdup.grouped.target_coverage",
                ".markdup.grouped.hs_metrics", ".markdup.grouped.qual_score_dist"]
    elif step == 2:
        return [".cram", ".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.stats", ".g.vcf.gz"]
    elif step == 3:
        return [".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.filteringStats.tsv"]


def submit_workflow(config, db, queue, workflow, sample=None, cohort=None, samples=None):
    '''
    Submits a workflow for execution

    COMMON: Refactorizar en una unica funcion la lectura de los JSON y envio a Cromwell

            workflow_json y genome_json son similares. Refactorizar a una unica funcion

    COMMON: Enviar datos del Backend que se usara en Cromwell (igual que accounting o queue).

    COMMON: Simplificar el jSON inputs que se envia a Cromwell

            Algunos datos estan estrechament relacionados (p.e: todos los asociados con Runtime o con los directorios
            RAW, REFERENCES, etc... Podríamos encapsularlos para que el objeto que se submit sea mas simple y tb.
            simplificaria el paso de los parametros dentro del WDL

    :param runid:
    :param workflow:
    :param sample:
    :param queue:
    :return:
    '''

    def wdl():
        # TODO: Usar bloque try-catch en lugar de if-the-else (como workflow_json)
        wdl_file = Path(f'./wdl/{workflow}.wdl')

        if wdl_file.exists():
            with open(wdl_file, "rb") as h_wdl_file:
                result = dict(wdl_file=BytesIO(h_wdl_file.read()))
                ''' TODO: Lo hago de esta forma tan barroca (en lugar de pasar el path porque en un futuro me 
                    gustaria crear el wdl on-the-fly'''
        else:
            print(f"ERROR: can't find {workflow}")
            exit(1)

        return result

    def options():
        '''
            TODO: Definir <outdir>, <calldior>, <logdir>, etc.., a nivel global de al App (objeto Execution o similar)
        '''
        if workflow == 'cohort_candidates':
            rootdir = f"{config.get('DIR', 'analysis_dir')}/{cohort}"
        else:
            rootdir = chip_cfg.root_dir(config, sample=sample, workflow=workflow)
        return dict(options_file=BytesIO(json.dumps({
            "use_relative_output_paths": "true",
            "final_workflow_outputs_dir": f"{rootdir}/_out",
            "final_workflow_log_dir": f"{rootdir}/_log",
            "final_call_logs_dir": f"{rootdir}/_call"
        }).encode()))

    def inputs():

        def workflow_json():
            '''
            :return: JSON corresponding to workflow submited
            '''

            json_file = Path(f'./wdl/{workflow}.json')

            try:
                with open(json_file, "rb") as h_json_file:
                    json_text = h_json_file.read().decode("utf-8")
                    result = re.sub(r"[\n\t]*", "", json_text)
            except FileNotFoundError:
                print(f"ERROR: can't find JSON for {workflow}")
                exit(1)

            return result

        def genome_json():
            '''
            :return: JSON corresponding to the reference genome
            '''

            genome_metadata = config.get('DIR', 'genome_metadata')
            try:
                with open(genome_metadata, "rb") as h_json_file:
                    json_text = h_json_file.read().decode("utf-8")
                    result = re.sub(r"[\n\t]*", "", json_text)
            except FileNotFoundError:
                print(f"ERROR: can't find genome metadata")
                exit(1)

            return result

        def create_inputs_json():
            """
            :return: JSON for inputs corresponding to the workflow

            NOTA:
                - OJO. La ruta de germline_calling esta cableada
            """

            def consensus_reads():
                def reads_in():
                    """ Returns a list with path to sample uBAM files """

                    def get_uBAM_path(json_str):
                        """ Generates de path of the uBAM file """

                        data = json.loads(json_str)
                        return f"""{config.get('DIR', 'raw_dir')}/{data.get('uBAM_path')}/{data.get('uBAM')}"""

                    try:
                        # db = chip_cfg.connect_db(config)
                        query = f"""SELECT json FROM samples WHERE sample = '{sample}' AND inclusion='enabled'"""
                        # db.cursor().execute("SELECT json FROM samples WHERE sample = '?' AND runid = '?'", sample, runid)
                        # result = json.JSONDecoder().decode("".join(cur.fetchone()))
                        # result_tmp = '","'.join([get_uBAM_path(json_item[0]) for json_item in chip_cfg.exec_query(db, query)])
                        result_tmp = '","'.join(
                            [get_uBAM_path(json_item[0]) for json_item in chip_cfg.exec_query(db, query)])
                        return f'["{result_tmp}"]'
                    except Exception as e:
                        print(f"ERROR: can't get {sample} info")
                        exit(1)
                    return result

                json_dict = f'''{{
                    "{workflow}.sample_id": "{sample}",
                    "{workflow}.rootdir": "{chip_cfg.root_dir(config, workflow=workflow, sample=sample)}",
                    "{workflow}.queue": "{queue}",
                    "{workflow}.accounting": "{config.get('ACCOUNTING', 'accounting')}",
                    "{workflow}.wf": {workflow_json()},
                    "{workflow}.genome": {genome_json()},
                    "{workflow}.raw_dir": "{config.get('DIR', 'raw_dir')}",
                    "{workflow}.references_dir": "{config.get('DIR', 'references_dir')}",
                    "{workflow}.bin_dir": "{config.get('DIR', 'bin_dir')}",
                    "{workflow}.enrichment_intervals": "{config.get('DIR', 'enrichment_intervals')}",
                    "{workflow}.reads_in": {reads_in()}
                }}'''

                return json_dict.encode()

            def variant_discovery():
                json_dict = f'''{{
                    "{workflow}.sample_id": "{sample}",
                    "{workflow}.queue": "{queue}",
                    "{workflow}.accounting": "{config.get('ACCOUNTING', 'accounting')}",
                    "{workflow}.wf": {workflow_json()},
                    "{workflow}.genome": {genome_json()},
                    "{workflow}.references_dir": "{config.get('DIR', 'references_dir')}",
                    "{workflow}.intervals_fullpath": "{config.get('DIR', 'enrichment_intervals')}",
                    "{workflow}.gnomad": "gnomAD/af-only-gnomad.hg38.vcf.gz",
                    "{workflow}.make_bamout": "False",
                    "{workflow}.run_ob_filter": "False",
                    "{workflow}.alignment_in": "{chip_cfg.root_dir(config, workflow='consensus_reads', sample=sample)}/_out/{sample}.bam"
                }}'''

                return json_dict.encode()

            def somatic_calling():
                json_dict = f'''{{
                    "{workflow}.sample_id": "{sample}",
                    "{workflow}.queue": "{queue}",
                    "{workflow}.accounting": "{config.get('ACCOUNTING', 'accounting')}",
                    "{workflow}.wf": {workflow_json()},
                    "{workflow}.genome": {genome_json()},
                    "{workflow}.references_dir": "{config.get('DIR', 'references_dir')}",
                    "{workflow}.vcf_in": "{chip_cfg.root_dir(config, workflow='variant_discovery', sample=sample)}/_out/{sample}.vcf.gz"
                }}'''

                return json_dict.encode()

            def cohort_candidates():
                json_dict = f'''{{
                    "{workflow}.cohort_id": "{cohort}",
                    "{workflow}.samples": {samples},
                    "{workflow}.queue": "{queue}",
                    "{workflow}.accounting": "{config.get('ACCOUNTING', 'accounting')}",
                    "{workflow}.wf": {workflow_json()},
                    "{workflow}.genome": {genome_json()},
                    "{workflow}.references_dir": "{config.get('DIR', 'references_dir')}"
                }}'''

                return json_dict.encode()

            def germline_calling():
                json_dict = f'''{{
                    "{workflow}.sample_id": "{sample}",
                    "{workflow}.rootdir": "{chip_cfg.root_dir(config, workflow=workflow, sample=sample)}",
                    "{workflow}.queue": "{queue}",
                    "{workflow}.accounting": "{config.get('ACCOUNTING', 'accounting')}",
                    "{workflow}.wf": {workflow_json()},
                    "{workflow}.genome": {genome_json()},
                    "{workflow}.references_dir": "{config.get('DIR', 'references_dir')}",
                    "{workflow}.intervals_fullpath": "{config.get('DIR', 'enrichment_intervals')}",
                    "{workflow}.gnomad": "gnomAD/af-only-gnomad.hg38.vcf.gz",
                    "{workflow}.alignment_in": "/data3/200831_JFuster_Somatic_PESA/_consensus/{sample}.bam"
                }}'''
                return json_dict.encode()

            switcher = {
                'consensus_reads': consensus_reads,
                'variant_discovery': variant_discovery,
                'somatic_calling': somatic_calling,
                'cohort_candidates': cohort_candidates,
                'germline_calling': germline_calling
            }

            inputs_func = switcher.get(workflow)

            return inputs_func()

        inputs = create_inputs_json()
        return dict(inputs_files=BytesIO(inputs))

    cromwell_server = config.get('CROMWELL', 'cromwell_server')
    stderr = StringIO()  # To avoid waring message for nor authenticating
    with contextlib.redirect_stderr(stderr):
        auth = CromwellAuth.harmonize_credentials(url=cromwell_server)
    command = getattr(CromwellAPI, 'submit')
    with contextlib.redirect_stderr(stderr):
        API_result = command(auth, **wdl(), **options(), **inputs())
    if isinstance(API_result, requests.Response):
        result = json.loads(API_result.text)
        if result['status'] == 'Submitted':
            return json.loads(API_result.text)['id']
        else:
            print(f"ERROR: workflow can't be submitted - It is in state {result['status']}")
            return None
    else:
        print(f'ERROR: {API_result}')
        exit(1)


def track(config, db, state, sample=None, step=None, uuid=None):
    """ Tracks sample executions

        TODO: Prodria modificarse para que no se requiera pasar el step_id ; porque asumimos que esta funcion
             no elimina trazas uy por tanto Tracks deberia tener el ultimo paso. Para eliminar trazas podría
             implementarse una funcion adicionañ
    """

    project = config.get('ACCOUNTING', 'project_name')

    if sample and step:
        query = f"REPLACE INTO Tracker(project, sample, step_idx, uuid, state) " \
                f"VALUES('{project}', '{sample}', '{step}', '{uuid}', '{state}')"
    elif not sample and not step and uuid:
        query = f"UPDATE Tracker SET state='{state}' WHERE uuid='{uuid}'"
    else:
        print(f"ERROR: bad call to 'track' function")
        exit(1)
    try:
        if is_valid_state(state):
            chip_cfg.exec_query(db, query)
    except Exception as e:
        print(f"ERROR: can't access track database : {str(e)}")


# Handlers
def init(config):
    """ Initialize DB tables

        TODO: Añadir indices a la BD
    """

    try:
        db = config.get('DB', 'db')
        with sqlite3.connect(db) as conn:

            # Samples
            conn.execute('DROP TABLE IF EXISTS samples')
            conn.execute("""CREATE TABLE samples(
                sample TEXT,
                runid TEXT,
                json TEXT,
                inclusion TEXT,
                PRIMARY KEY (sample, runid))""")

            # Tracker
            conn.execute('DROP TABLE IF EXISTS Tracker')
            conn.execute("""CREATE TABLE Tracker(
                sample TEXT,
                step_idx INTEGER,
                uuid TEXT,
                state TEXT,
                PRIMARY KEY (sample, step_idx))""")

            """
                OTA: Para facilitar la implementancion de la vista Tracks ; se establece como PK (sample, step_idx) 
                de modo que una nueva ejecucion sobreescribe la info (uuid, state) de la ejecucion anterior.
            """

            # Tracks (Current track of a sample)
            conn.execute('DROP VIEW IF EXISTS Tracks')
            # SQLITE
            # conn.execute("""CREATE VIEW Tracks AS
            #     SELECT T.sample, T.step_idx as current_step, T.uuid, T.state FROM Tracker AS T
            #     WHERE T.sample IN (SELECT sample FROM samples WHERE inclusion='enabled')
            #     GROUP BY T.sample
            #     HAVING max(T.step_idx)""")
            conn.execute("""CREATE VIEW CHIP_INSPECTOR.Tracks AS SELECT project, sample, max_step_idx as current_step
             FROM CHIP_INSPECTOR.max_step AS M
             WHERE M.sample in (SELECT sample FROM CHIP_INSPECTOR.samples WHERE inclusion='enabled'""")

            # La misma vista para MYSQL
            mysql = """CREATE OR REPLACE VIEW CHIP_INSPECTOR.max_step AS
                SELECT project, sample, MAX(step_idx) as max_step_idx FROM CHIP_INSPECTOR.Tracker GROUP BY project, sample"""
            mysql = """CREATE OR REPLACE VIEW CHIP_INSPECTOR.Tracks AS
                        SELECT t1.project, t1.sample, t1.step_idx as current_step, t1.uuid, t1.state FROM CHIP_INSPECTOR.Tracker t1
                        INNER JOIN CHIP_INSPECTOR.max_step t2 ON t1.project = t2.project AND t1.sample = t2.sample AND t1.step_idx = t2.max_step_idx"""

            # Cohort
            conn.execute('DROP TABLE IF EXISTS Cohort')
            conn.execute("""CREATE TABLE Cohort(
                cohort TEXT,
                uuid TEXT,
                state TEXT,
                PRIMARY KEY (cohort))""")

            conn.commit()

            # Create tables por 'base' stats
            db_stats = config.get('DB', 'db_stats')
            stats.create_base_stats(db_stats)

    except Exception as e:
        print(f"ERROR: can't initialize project : {str(e)}")
        exit(1)


def load(config, db, samplesheet_file, samplesheet_project, runid, flag_skip_missing_samples=False):
    """ Loads samples into a project"""

    def to_json():
        """ Generate JSON from SampleSheet. """

        def samplesheet():
            """ Returns a dict sample : info ; where info is a list dict (one for each row)."""

            def sample_name():
                """ Get sample name from SampleSheet field.

                    TODO: Permitir la obtencion del nombre de muestra sea configurable por el usuario
                        Por ejemplo mediante regexp o similar
                """
                return line['Sample_Name'].split('_', 1)[1]
                # return line['Sample_Name']

            def parse():
                """ Returns a list of dictionaries with the sample info parsed from 'samplesheet'. """
                result = []
                with open(samplesheet_file) as sheet:
                    flag = False
                    for row in csv.reader(sheet):
                        if row[0] == '[Data]':  # A partir de aqui empezamos a leer las filas de los ficheros fastq
                            header = next(csv.reader(sheet))
                            flag = True
                        elif flag:
                            result.append(dict(zip(header, row)))
                return result

            result = dict()
            for line in parse():
                if line.get('Sample_Project', None) == samplesheet_project:
                    result.setdefault(sample_name(), []).append(line)

            return result

        # Transform to JSON-like from Samplesheet-like
        def json_like(sample, sample_info):
            def mutate_2_json(file_info):

                # Returns path
                def path():
                    def run_folder():
                        return f"{runid}"

                    def project_folder():
                        return f"{file_info['Sample_Project']}"

                    def sample_folder():
                        return f"{file_info['Sample_Name']}"

                    return f"{run_folder()}/{project_folder()}/{sample_folder()}"

                # Return uBAM
                #   Requires just one uBam per Samplesheet row
                def uBam():
                    pattern = f"{file_info['Sample_Name']}*_unmapped.bam"
                    try:
                        return next(Path(f"{config.get('DIR', 'raw_dir')}/{path()}").glob(pattern)).name

                    except StopIteration:
                        return None

                # Return preparation pool
                def pool():
                    return file_info["Sample_ID"].split('_', 2)[0]

                # Return preparation well
                def well():
                    return file_info["Sample_ID"].split('_', 2)[1]

                return {
                    "name": f"{sample}",
                    "samplesheet_project": f"{file_info['Sample_Project']}",
                    "uBAM_path": path(),
                    "uBAM": uBam(),
                    "pool": pool(),
                    "well": well()
                }

            return mutate_2_json(sample_info)

        return {sample: json_like(sample, info[0]) for (sample, info) in samplesheet().items()}

    def write_json():
        """
            Write JSON to DB.
        """

        def add_row():
            try:
                project = config.get('ACCOUNTING', 'project_name')
                sql = "REPLACE INTO samples(project, sample, runid, json, inclusion, experiment) " \
                      "VALUES('{}', '{}', '{}', '{}', 'enabled', '{}')" \
                    .format(project, sample, runid, json.dumps(json_text), json_text.get('samplesheet_project'))
                db.cursor().execute(sql)
            except Exception as e:
                raise RuntimeError(f"ERROR: can't add data from {runid} of {sample} : {str(e)}")

        add_row()

    for sample, json_text in to_json().items():
        if json_text.get('uBAM', None) is not None:
            write_json()
            # TODO Comprobar si la carga ha sido existosa hacer el track correspoendiente
            # Debemos asegurarnos que sample y runid no son null
            # Asumimos, por el momento, que es correcta y por tanto state = '_ready'
            track(config, db, '_ready', step=1, sample=sample)
        else:
            if flag_skip_missing_samples:
                print(f"WARNING: No uBAM for {sample}")
            else:
                print(f"ERROR: No uBAM for {sample}")
                exit(1)


def BAK_load(config, db, samplesheet_file, samplesheet_project, runid, flag_skip_missing_samples=False):
    """ Loads samples into a project"""

    def to_json():
        """ Generate JSON from SampleSheet. """

        def samplesheet():
            """ Returns a dict sample : info ; where info is a list dict (one for each row)."""

            def sample_name():
                """ Get sample name from SampleSheet field.

                    TODO: Permitir la obtencion del nombre de muestra sea configurable por el usuario
                        Por ejemplo mediante regexp o similar
                """
                return line['Sample_Name'].split('_', 1)[1]
                # return line['Sample_Name']

            def parse():
                """ Returns a list of dictionaries with the sample info parsed from 'samplesheet'. """
                result = []
                with open(samplesheet_file) as sheet:
                    flag = False
                    for row in csv.reader(sheet):
                        if row[0] == '[Data]':  # A partir de aqui empezamos a leer las filas de los ficheros fastq
                            header = next(csv.reader(sheet))
                            flag = True
                        elif flag:
                            result.append(dict(zip(header, row)))
                return result

            result = dict()
            for line in parse():
                if line.get('Sample_Project', None) == samplesheet_project:
                    result.setdefault(sample_name(), []).append(line)

            return result

        # Transform to JSON-like from Samplesheet-like
        def json_like(sample, sample_info):
            def mutate_2_json(file_info):

                # Returns path
                def path():
                    def run_folder():
                        return f"{runid}"

                    def project_folder():
                        return f"{file_info['Sample_Project']}"

                    def sample_folder():
                        return f"{file_info['Sample_Name']}"

                    return f"{run_folder()}/{project_folder()}/{sample_folder()}"

                # Return uBAM
                #   Requires just one uBam per Samplesheet row
                def uBam():
                    pattern = f"{file_info['Sample_Name']}*_unmapped.bam"
                    try:
                        return next(Path(f"{config.get('DIR', 'raw_dir')}/{path()}").glob(pattern)).name

                    except StopIteration:
                        if flag_skip_missing_samples:
                            print(f"WARNING: No fastq for {path()}/{pattern}")
                            return None
                        else:
                            print(f"ERROR: No fastq for {path()}/{pattern}")
                            exit(1)

                # Return preparation pool
                def pool():
                    return file_info["Sample_ID"].split('_', 2)[0]

                # Return preparation well
                def well():
                    return file_info["Sample_ID"].split('_', 2)[1]

                return {
                    "name": f"{sample}",
                    "uBAM_path": path(),
                    "uBAM": uBam(),
                    "pool": pool(),
                    "well": well()
                }

            return mutate_2_json(sample_info)

        return {sample: json.dumps(json_like(sample, info[0])) for (sample, info) in samplesheet().items()}

    def write_json():
        """
            Write JSON to DB.
        """

        def add_row():
            try:
                project = config.get('ACCOUNTING', 'project_name')
                sql = "REPLACE INTO samples(project, sample, runid, json, inclusion) VALUES('{}', '{}', '{}', '{}', 'enabled')" \
                    .format(project, sample, runid, json_text)
                db.cursor().execute(sql)
            except Exception as e:
                raise RuntimeError(f"ERROR: can't add data from {runid} of {sample} : {str(e)}")

        add_row()

    for sample, json_text in to_json().items():
        write_json()
        # TODO Comprobar si la carga ha sido existosa hacer el track correspoendiente
        # Debemos asegurarnos que sample y runid no son null
        # Asumimos, por el momento, que es correcta y por tanto state = '_ready'
        track(config, db, '_ready', step=1, sample=sample)


def drop(config, runid, sample):
    try:
        db = config('DB', 'db')
        with sqlite3.connect(db) as conn:

            # Samples
            conn.execute(f"UPDATE FROM 'samples' "
                         f"SET inclusion = 'dropped' "
                         f"WHERE sample = '{sample}' AND runid = '{runid}'")
            conn.commit()
            conn.close()
    except Exception as e:
        print(f"ERROR: can't access sample data : {str(e)}")


def enable(config, runid, sample):
    try:
        db = config('DB', 'db')
        with sqlite3.connect(db) as conn:

            # Samples
            conn.execute(f"UPDATE FROM samples "
                         f"SET enabled = 'enabled' "
                         f"WHERE sample = '{sample}' AND runid = '{runid}'")
            conn.commit()
            conn.close()
    except Exception as e:
        print(f"ERROR: can't access sample data : {str(e)}")


def disable(config, runid, sample):
    try:
        db = config('DB', 'db')
        with sqlite3.connect(db) as conn:

            # Samples
            conn.execute(f"UPDATE FROM 'samples' "
                         f"SET enabled = 'disabled' "
                         f"WHERE sample = '{sample}' AND runid = '{runid}'")
            conn.commit()
            conn.close()
    except Exception as e:
        print(f"ERROR: can't access sample data : {str(e)}")


def create_panel(config, panel_name, panel_file):
    try:
        # TODO: Leer Conector del BD del .ini
        db_connection = TargetPanelMySQLConnector(config)
        TargetPanel.create_from_file(config, db_connection, panel_name, panel_file)
    except Exception as e:
        print(f"ERROR: can't create panel : {str(e)}")


def remove_panel(config, panel_name):
    try:
        # TODO: Leer Conector del BD del .ini
        db_connection = TargetPanelMySQLConnector(config)
        TargetPanel(config, db_connection, panel_name).delete()
    except Exception as e:
        print(f"ERROR: can't remove panel : {str(e)}")


def process(config, db, sample, workflow, flag_test):
    queue = config.get('ACCOUNTING', 'test_queue') if flag_test else config.get('ACCOUNTING', 'working_queue')
    uuid = submit_workflow(config, db, queue=queue, workflow=workflow[0], sample=sample)
    if uuid:
        print(f"Running {workflow[0]} for {sample} as job {uuid}")
        track(config, db, state='_in_progress', sample=sample, step=workflow[1], uuid=uuid)
    db.commit()


def abort(config, db, uuid):
    stderr = StringIO()  # To avoid waring message for nor authenticating
    with contextlib.redirect_stderr(stderr):
        cromwell_server = config.get('CROMWELL', 'cromwell_server')
        auth = CromwellAuth.harmonize_credentials(url=cromwell_server)
    command = getattr(CromwellAPI, 'abort')
    with contextlib.redirect_stderr(stderr):
        API_result = command(auth=auth, uuid=uuid)
        # No compruebo salida porque si el uuid no esta devuelve 404. Habria que hacer algo mas depurado
        if isinstance(API_result, requests.Response):
            result = json.loads(API_result.text)
            if not result['status'] == 'Aborting':
                print(f"WARNING: Cromwell have not aborted  {uuid}."
                      f"However CHIP-INSPECTOR marked as aborted.")
            print(f"Aborting {uuid} results '{result['status']}'")
            track(config, db, state='_aborted', uuid=uuid)
        else:
            print("ERROR")
            exit(1)


def check(config, db, project_name, sample, current_step, uuid):
    def cromwell_workflow_status():
        """ Get Workflow state from cromwell DB

        """
        try:
            cromwell_host = config.get('CROMWELL', 'cromwell_host')
            cromwell_port = int(config.get('CROMWELL', 'cromwell_port'))
            cromwell_user = config.get('CROMWELL', 'cromwell_user')
            cromwell_passwd = config.get('CROMWELL', 'cromwell_passwd')
            cromwell_db = config.get('CROMWELL', 'cromwell_db')
            conn = pymysql.connect(host=cromwell_host, port=cromwell_port, user=cromwell_user, passwd=cromwell_passwd,
                                   db=cromwell_db)

            with conn.cursor() as cur:
                sql = "SELECT METADATA_VALUE FROM cromwell.METADATA_ENTRY " \
                      "WHERE WORKFLOW_EXECUTION_UUID=%s and METADATA_KEY='status' " \
                      "and METADATA_TIMESTAMP=(SELECT max(METADATA_TIMESTAMP) from cromwell.METADATA_ENTRY " \
                      "WHERE WORKFLOW_EXECUTION_UUID=%s " \
                      "and METADATA_KEY='status')"
                cur.execute(sql, (uuid, uuid))
                result = cur.fetchone()
        except Exception as e:
            print(f"ERROR: can't check state of workflow {uuid} : {str(e)}")
            exit(1)
        finally:
            conn.close()

        if result is None:
            print(f"ERROR: No workflow {uuid} accounted in Cromwell")
        else:
            return result[0]

    def delete_cromwell_execution_dir():
        dir = "{}/{}/{}".format(config.get('CROMWELL', 'cromwell_executions'), get_workflow(current_step), uuid)
        shutil.rmtree(dir, ignore_errors=True, onerror=None)

    def load_basic_cov_stats(config, db, project_name, sample):
        def add_sample_cov_stat(project_name, sample, df_raw, df_consensus):
            """ Añade info sobre cobertura media global exceptuando chr y chrX para evitar el sesgo por genero """
            # No tenemos en cuenta los chromosomas X y Y para que no se vea sesgado por el sexo del indiviuo
            avg_cov_raw = df_raw.query("chrom not in ('chrX','chrY')")['coverage'].mean()
            avg_cov_consensus = df_consensus.query("chrom not in ('chrX','chrY')")['coverage'].mean()
            ratio_cov = avg_cov_consensus / avg_cov_raw if avg_cov_raw > 0 else 0

            query = f"REPLACE INTO CHIP_INSPECTOR.sample_coverage(project, sample, cov_raw, cov_consensus, ratio_cov) " \
                    f"VALUES ('{project_name}', '{sample}', {avg_cov_raw}, {avg_cov_consensus}, {ratio_cov})"
            try:
                chip_cfg.exec_query(db, query)
            except Exception as e:
                print(f"ERROR: adding data to database for {project_name} - {sample}: {str(e)}")

        def add_sample_cov_stat_by_gene(project_name, sample, df_raw, df_consensus):
            """ Añade info con cobertura por gen"""
            avg_cov_raw_by_gene = df_raw.groupby(['target'])['coverage'].mean().reset_index(['target'])
            avg_cov_consensus_by_gene = df_consensus.groupby(['target'])['coverage'].mean().reset_index(['target'])
            df_sample = pd.merge(avg_cov_raw_by_gene, avg_cov_consensus_by_gene, how='outer', on='target',
                                 suffixes=('_raw', '_consensus'))
            df_sample['sample'] = sample

            try:
                for i, row in df_sample.iterrows():
                    cov_raw = float(row['coverage_raw'])
                    cov_consensus = float(row['coverage_consensus'])
                    gene = row['target']
                    ratio_cov = cov_consensus / cov_raw if cov_raw > 0 else 0
                    query = f"REPLACE INTO CHIP_INSPECTOR.sample_coverage_by_gene(project, sample, gene, cov_raw, " \
                            f"cov_consensus, ratio_cov) VALUES ('{project_name}', '{sample}', '{gene}', {cov_raw}, " \
                            f"{cov_consensus}, {ratio_cov})"
                    chip_cfg.exec_query(db, query)
            except Exception as e:
                print(f"ERROR: adding data to database for {project_name} - {sample} : {str(e)}")

        basedir = f"{config.get('DIR', 'analysis_dir')}/{sample}/consensus_reads/_out/"
        df_raw = pd.read_csv(f'{basedir}/{sample}.markdup.grouped.base_coverage', sep='\t', header=0)
        df_consensus = pd.read_csv(f'{basedir}/{sample}.base_coverage', sep='\t', header=0)

        add_sample_cov_stat_by_gene(project_name, sample, df_raw, df_consensus)
        add_sample_cov_stat(project_name, sample, df_raw, df_consensus)
        db.commit()

    def succeeded():
        """ Checks if succeeded

            COMMON: Adds post-condition checking
        """

        def exist_result_files(suffixes):
            """ Checks if all result files exists"""

            def exist_file(suffix):
                """ Returns the path for the file"""
                outdir = f"{chip_cfg.root_dir(config, workflow=workflow, sample=sample)}/_out"
                return Path(f"{outdir}/{sample}{suffix}").exists()

            return reduce(lambda x, y: x and y, [exist_file(suffix) for suffix in suffixes])


        suffixes = get_result_files_suffixes(current_step)
        # if exist_result_files(suffixes) and stats.populate_stats(workflow, sample, runid):
        if exist_result_files(suffixes):
            # BORRAMOS la ejecucion de cromwell-executions
            delete_cromwell_execution_dir()
            if current_step == 1:
                load_basic_cov_stats(config, db, project_name, sample)
            return ('_successful', None)
        else:
            return ('_failed', 'result files do not exist')

    def running():
        return ("_in_progress", None)

    def error():
        return ('_failed', 'Workflow status in cromwell labeled as error')

    def aborted():
        return ('_aborted', None)

    def unknown():
        return ('_failed', 'Workflow status in cromwell not recognized')

    switcher = {
        'Succeeded': succeeded,
        'Submitted': running,
        'Running': running,
        'Failed': error,
        'Aborting': aborted,
        'Aborted': aborted
    }


    workflow = get_workflow(current_step)
    func = switcher.get(cromwell_workflow_status(), unknown)
    state, msg = func()

    if is_valid_state(state):
        if state == '_failed':
            print(f"Checking {sample} results '{state}' for {workflow}: {msg}")
        else:
            print(f"Checking {sample} results '{state}' for {workflow}")
        track(config, db, state=state, sample=sample, step=current_step, uuid=uuid)
        db.commit()
    else:
        print(f"ERROR: Checking {sample} results a unidentified '{state}'")


def cohort(db, project, cohort, sample):
    sql = f"REPLACE INTO cohorts(project, cohort, sample) " \
          f"VALUES ('{project}', '{cohort}', '{sample}')"
    db.cursor().execute(sql)


def report(config):
    """ Displays results.

    NOTA: Implementar un parametro --add para que el report 'somatic_call' se añada (y no substituya al anterior).


    """
    # rp.display()


def candidates(config, db, cohort, samples, flag_test):
    """
        Para obtener el listado de mutaciones candidatas se realian dos pasos:
        1) Se lanza el WDL 'cohort_candidates' que combina ls VCF de todas las muestas de la cohorte, normaliza el VCF
        y lo anota.
        2) Se parsea el VCF anotado para generar el fichero de mutaciones candidatas
    """
    queue = config.get('ACCOUNTING', 'test_queue') if flag_test else config.get('ACCOUNTING', 'working_queue')
    uuid = submit_workflow(config, db, queue=queue, workflow='cohort_candidates', cohort=cohort, samples=samples)
    if uuid:
        print(f"Preparing cohort {cohort} as job {uuid}")
        # Habria que tracker su ejecucion (tabla 'cohort', ver init)

    # TODO: Implementar la comprobacion del resultado


def fingerprint(config, sample, flag_test):
    queue = config.get('ACCOUNTING', 'test_queue') if flag_test else config.get('ACCOUNTING', 'working_queue')
    uuid = submit_workflow(config, queue=queue, workflow='germline_calling', sample=sample)
    if uuid:
        print(f"Calling germinal variants for {sample}")
        # Habria que tracker su ejecucion

def stat_cov(config, db, cohort, stat_outfile):
    """ Calcula estadisticas de cobertura para las muestras de una cohorte """
    def get_sample_cov_stats():
        panel = config.get('GENERAL', 'target_panel')
        query = f"SELECT sample, gene, cov_consensus FROM cohort_members cm " \
                f"INNER JOIN sample_coverage_by_gene scbg " \
                f"ON cm.experiment_id = scbg.project AND cm.sample_id = scbg.sample " \
                f"WHERE cm.cohort_id = '{cohort}' "\
                f"AND gene IN (SELECT gene FROM CHIP_INSPECTOR.target_panel tp WHERE panel_name = '{panel}')"
        try:
            #return chip_cfg.exec_query(db, query)
            return pd.DataFrame(chip_cfg.exec_query(db, query), columns=['sample', 'gene', 'cov'])
        except Exception as e:
            print(f"ERROR: can't obtain sample coverage data : {str(e)}")
            exit(1)

    df = get_sample_cov_stats()
    k1 = df.groupby(['sample'])['cov'].mean().reset_index(['sample'])
    k2 = k1['cov'].mean()
    pass
def dummy_load_basestats(config, db, project_name, sample, conversion_table):
    def add_sample_cov_stat(project_name, sample, df_raw, df_consensus):
        """ Añade info sobre cobertura media global exceptuando chr y chrX para evitar el sesgo por genero """
        # No tenemos en cuenta los chromosomas X y Y para que no se vea sesgado por el sexo del indiviuo
        avg_cov_raw = df_raw.query("chrom not in ('chrX','chrY')")['coverage'].mean()
        avg_cov_consensus = df_consensus.query("chrom not in ('chrX','chrY')")['coverage'].mean()
        ratio_cov = avg_cov_consensus / avg_cov_raw if avg_cov_raw > 0 else 0

        query = f"REPLACE INTO CHIP_INSPECTOR.sample_coverage(project, sample, cov_raw, cov_consensus, ratio_cov) " \
                f"VALUES ('{project_name}', '{sample}', {avg_cov_raw}, {avg_cov_consensus}, {ratio_cov})"
        try:
            chip_cfg.exec_query(db, query)
        except Exception as e:
            print(f"ERROR: adding data to database for {project_name} - {sample}: {str(e)}")

    def add_sample_cov_stat_by_gene(project_name, sample, df_raw, df_consensus):
        """ Añade info con cobertura por gen"""
        avg_cov_raw_by_gene = df_raw.groupby(['target'])['coverage'].mean().reset_index(['target'])
        avg_cov_consensus_by_gene = df_consensus.groupby(['target'])['coverage'].mean().reset_index(['target'])
        df_sample = pd.merge(avg_cov_raw_by_gene, avg_cov_consensus_by_gene, how='outer', on='target',
                             suffixes=('_raw', '_consensus'))
        df_sample['sample'] = sample

        try:
            for i, row in df_sample.iterrows():
                cov_raw = float(row['coverage_raw'])
                cov_consensus = float(row['coverage_consensus'])
                gene = row['target']
                ratio_cov = cov_consensus / cov_raw if cov_raw > 0 else 0
                query = f"REPLACE INTO CHIP_INSPECTOR.sample_coverage_by_gene(project, sample, gene, cov_raw, " \
                        f"cov_consensus, ratio_cov) VALUES ('{project_name}', '{sample}', '{gene}', {cov_raw}, " \
                        f"{cov_consensus}, {ratio_cov})"
                chip_cfg.exec_query(db, query)
        except Exception as e:
            print(f"ERROR: adding data to database for {project_name} - {sample} : {str(e)}")

    basedir = f"{config.get('DIR', 'analysis_dir')}/{sample}/consensus_reads/_out/"
    df_raw = pd.read_csv(f'{basedir}/{sample}.markdup.grouped.base_coverage', sep='\t', header=0)
    df_consensus = pd.read_csv(f'{basedir}/{sample}.base_coverage', sep='\t', header=0)

    # Convertimos el ID de la muestra para corregir muestras cambiadas (mal etiquetadas)
    new_sample = conversion_table.get(sample, sample)

    add_sample_cov_stat_by_gene(project_name, new_sample, df_raw, df_consensus)
    add_sample_cov_stat(project_name, new_sample, df_raw, df_consensus)
    db.commit()

#######################3

# def qc(workflow=main_workflow):
#     '''
#         COMMON: Trabajar a nivel de dataset (flowcell o user-defined) en lugar de workflow
#     '''
#     def samples():
#         def deserialize(name, stat, sdata):
#             # name, stat, sdata = row
#             return name, json.loads(stat), json.loads(sdata)
#
#         query = "SELECT T.sample as name, S.json as stat , P.json AS sample " \
#                 "FROM Tracker as T, 'Stats' as S  , 'samples' AS P " \
#                 "WHERE T.workflow=? and T.state='_successful' " \
#                     "AND S.workflow=T.workflow AND S.sample=T.sample " \
#                     "AND P.sample=T.sample"
#
#         with sqlite3.connect(db) as conn:
#             cur = conn.cursor()
#             cur.execute(query, (workflow,))
#             db_result = cur.fetchall()
#
#         if len(db_result) == 0:
#             print(f'ERROR: There is no samples of {workflow} for QC')
#             exit(1)
#         else:
#             return [deserialize(*row) for row in db_result]
#
#     def qc_table_reads(f):
#
#         '''
#             Para mirar la relacion entre las lecturas RAW y consensus (FilteredConsensus) por pool
#         '''
#
#
#         if header:
#             f.write(f"name\t"
#                     f"pool\t"
#                     f"raw\t"
#                     f"consensus\t"
#                     f"perc_dup\n"
#                     )
#
#         consensus = stat['filter_stderr']['filtered_emitted'] # ya es int
#         raw = int(stat['markdup_metrics']['READ_PAIRS_EXAMINED']) * 2 \
#             + int(stat['markdup_metrics']['UNPAIRED_READS_EXAMINED']) \
#             + int(stat['markdup_metrics']['UNMAPPED_READS'])
#         perc_dup = float(stat['markdup_metrics']['PERCENT_DUPLICATION'])
#
#         f.write(f"{name}\t"
#               f"{sdata['pool']}\t"
#               f"{raw}\t"
#               f"{consensus}\t"
#               f"{perc_dup}\n"
#               )
#
#     def qc_table_umi(f):
#         '''
#             Para mirar el comportamiento de los UMI (umi_metrics)
#         '''
#
#         '''
#         if header:
#             print("".join("{}\t".format(k) for k in stat['umi_metrics'].keys()))
#         else:
#             print("".join("{}\t".format(k) for k in stat['umi_metrics'].values()))
#         '''
#
#         if header:
#             f.write(f"name\t"
#                   f"mean_umi_length\t"
#                   f"duplicate_sets_ignoring_umi\t"
#                   f"duplicate_sets_with_umi\t"
#                   f"observed_unique_umis\t"
#                   f"inferred_unique_umis\t"
#                   f"observed_umi_entropy\t"
#                   f"inferred_umi_entropy\n")
#
#         mean_umi_length = int(stat['umi_metrics']['MEAN_UMI_LENGTH'])
#         duplicate_sets_ignoring_umi = int(stat['umi_metrics']['DUPLICATE_SETS_IGNORING_UMI'])
#         duplicate_sets_with_umi = int(stat['umi_metrics']['DUPLICATE_SETS_WITH_UMI'])
#         observed_unique_umis = int(stat['umi_metrics']['OBSERVED_UNIQUE_UMIS'])
#         inferred_unique_umis = int(stat['umi_metrics']['INFERRED_UNIQUE_UMIS'])
#         observed_umi_entropy = float(stat['umi_metrics']['OBSERVED_UMI_ENTROPY'])
#         inferred_umi_entropy = float(stat['umi_metrics']['INFERRED_UMI_ENTROPY'])
#
#
#         f.write(f"{name}\t"
#               f"{mean_umi_length}\t"
#               f"{duplicate_sets_ignoring_umi}\t"
#               f"{duplicate_sets_with_umi}\t"
#               f"{observed_unique_umis}\t"
#               f"{inferred_unique_umis}\t"
#               f"{observed_umi_entropy}\t"
#               f"{inferred_umi_entropy}\n")
#
#     def qc_table_family_size(f):
#
#         if header:
#             f.write(f"name\t"
#                   f"family_size\t"
#                   f"fraction_gt_or_eq_family_size\n")
#
#         for bin in stat['grouped_hist']:
#             family_size = int(bin['family_size'])
#             fraction_gt_or_eq_family_size = float(bin['fraction_gt_or_eq_family_size'])
#
#             f.write(f"{name}\t"
#               f"{family_size}\t"
#               f"{fraction_gt_or_eq_family_size}\n")
#
#     def qc_table_qual_score_dist(f):
#
#         if header:
#
#             f.write(f"name\t"
#                   f"QUALITY\t"
#                   f"COUNT_OF_Q\n")
#
#         for bin in stat['qual_score_dist']:
#             quality = int(bin['QUALITY'])
#             count_of_q = int(bin['COUNT_OF_Q'])
#
#             f.write(f"{name}\t"
#                     f"{quality}\t"
#                     f"{count_of_q}\n")
#
#     def qc_table_coverage(f, pattern):
#                 '''
#                     Parte de esta funcion deberia implentarse en 'check' dentro de 'base_coverage'
#                 '''
#
#                 if header:
#                     f.write(f"name\t"
#                             f"chr\t"
#                             f"pos\t"
#                             f"target\t"
#                             f"cov\n")
#
#                 bc_path = next(Path(f'{root_dir(workflow, name)}/_out/').glob(pattern))
#                 with open(bc_path, 'r') as file_bc:
#                     flag = True
#                     for chrom, pos, target, coverage in csv.reader(file_bc, delimiter='\t'):
#                         if flag:
#                             flag = False
#                         else:
#                             f.write(f"{name}\t"
#                                     f"{chrom}\t"
#                                     f"{pos}\t"
#                                     f"{target}\t"
#                                     f"{coverage}\n")
#
#     base_dir = f'{analysis_dir}/{workflow}/_QC'
#     base_dir.mkdir(parents=True, exist_ok=True)
#
#     f1 = open(f'{base_dir}/qc_table_reads.tsv', 'w')
#     f2 = open(f'{base_dir}/qc_table_umi.tsv', 'w')
#     f3 = open(f'{base_dir}/qc_table_family_size.tsv', 'w')
#     f4 = open(f'{base_dir}/qc_table_qual_score_dist.tsv', 'w')
#     f5 = open(f'{base_dir}/qc_raw_coverage.tsv', 'w')   # TODO: Combinar f5 y f6
#     f6 = open(f'{base_dir}/qc_consensus_coverage.tsv', 'w')
#
#     header = True
#     for name, stat, sdata in samples():
#         qc_table_reads(f1)
#         qc_table_umi(f2)
#         qc_table_family_size(f3)
#         qc_table_qual_score_dist(f4)
#         qc_table_coverage(f5, '*.mapped.markdup.grouped.base_coverage')
#         qc_table_coverage(f6, '*.base_coverage')
#         header = False
#
#     f1.close()
#     f2.close()
#     f3.close()
#     f4.close()
#     f5.close()
#     f6.close()
#     pass
#
#
# def submit(workflow, runid, sample, flag_test):
#     uuid = submit_workflow(workflow=workflow, sample=sample, queue=f"{'test' if flag_test else 'prod'}")
#     track('_in_progress', mode='new_workflow', sample=sample, runid=runid, uuid=uuid, workflow=workflow)
