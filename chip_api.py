import csv
import os
import shutil
import sys

import chip_parsing_crossectional
import chip_parsing_longitudinal as chip_longitudinal
import chip_parsing_crossectional as chip_crossectional
from chip_cfg import *
from pathlib import Path
import chip_bk

# Es para test se puede quitar despues
import TargetPanel as tp


def final_process_step():
    """ Devuelve el ultimo 'step' del pipeline de procesado"""
    return ('somatic_calling', 3)


# Handlers
def init(config, flag_force):
    """ Initialize project.

    TODO: Implementar --force

        Actualmente siempre se pide --force pero en realidad solo deberia usarse y existe la bd y las tablas estan
        creadas; es decir, que el proyecto se hubiera inicializado anteriormente.

    """

    if flag_force:
        chip_bk.init(config)
    else:
        raise RuntimeError("You should use --force to overwrite data.")


def load(config, samplesheet, samplesheet_project, runid, flag_skip_missing_samples=False):
    """ Loads sequences of a sample from a flowcell

        COMMON: Manejar ambos uBAM y Samplesheet

        TODO: Inicializar una muesta/runid descartada con 'drop'
    """

    db = connect_db(config)

    if not Path(samplesheet).exists():
        raise FileNotFoundError(f"Samplesheet doesn't exist.")

    try:
        chip_bk.load(config, db, samplesheet, samplesheet_project, runid, flag_skip_missing_samples)
        db.commit()
    except Exception as e:
        db.rollback()
        print(f"ERROR: {str(e)}")
    finally:
        db.close()


def drop(config, runid, sample):
    """ Discards samples from the project.

        TODO: Implement 'drop'
    """
    pass


def enable(config, runid, sample):
    """ Enable samples.

        MAIN: Combinar en un unico comando enable/disable (en lugar de tener 2 comandos separados)

            Estos 2 comandos son casi iguales por lo que quiza sea buena idea combinarlos en uno solo.

    """

    def samples():
        """ Get samples to be processed"""

        if sample is not None and runid is not None:
            query = f"SELECT runid, sample  FROM 'Tracker' WHERE step_idx= 0 AND project = '{project_name} " \
                    f"AND sample = '{sample}' AND runid = '{runid}'"
        elif sample is not None:
            query = f"SELECT runid, sample  FROM 'Tracker' WHERE step_idx= 0 AND project = '{project_name} " \
                    f"AND sample = '{sample}'"
        elif runid is not None:
            query = f"SELECT runid, sample  FROM 'Tracker' WHERE step_idx= 0 AND project = '{project_name} " \
                    f"AND runid = '{runid}'"
        else:
            raise RuntimeError("No sample were specified.")

        try:
            return exec_query(db, query)
        except Exception as e:
            print(f"ERROR: can't read sample data : {str(e)}")

    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')
    for runid, sample in samples():
        chip_bk.enable(config, runid, sample)


def disable(config, runid, sample):
    """ Disable samples"""

    def samples():
        """ Get samples to be processed"""

        if sample is not None and runid is not None:
            query = f"SELECT runid, sample  FROM 'Tracker' WHERE step_idx= 0 AND project = '{project_name} " \
                    f"AND sample = '{sample}' AND runid = '{runid}'"
        elif sample is not None:
            query = f"SELECT runid, sample  FROM 'Tracker' WHERE step_idx= 0 AND project = '{project_name}" \
                    f" AND sample = '{sample}'"
        elif runid is not None:
            query = f"SELECT runid, sample  FROM 'Tracker' WHERE step_idx= 0 AND project = '{project_name} " \
                    f"AND runid = '{runid}'"
        else:
            raise RuntimeError("No sample were specified.")

        try:
            return exec_query(db, query)
        except Exception as e:
            print(f"ERROR: can't read sample data : {str(e)}")

    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')
    for runid, sample in samples():
        chip_bk.disable(config, runid, sample)


def create_panel(config, panel_name, panel_file):
    """ Creates a new target panel """
    if Path(panel_file).exists():
        chip_bk.create_panel(config, panel_name, panel_file)
    else:
        raise FileNotFoundError(f"Panel file doesn't exist.")


def remove_panel(config, panel_name):
    """ Deletes a target panel """
    chip_bk.remove_panel(config, panel_name)


def process(config, runid, sample, reset_idx, flag_test=False, batch_size=None):
    """  Process samples. """

    # TODO: Implementar '--reset'

    def delete_trace(config, sample, reset_idx):
        """ Delete traces of downstream analysis from r' set point """

        # TODO: Comprobar que las muestras no estan `in_progress`
        #       Aunque, en teoria, esta circunstancia no se da porque aqui solo llegan muestras
        #       que no estan en ese estado

        query = f"DELETE FROM Tracker WHERE project = '{project_name}' " \
                f"AND sample = '{sample}' AND step_idx >= {reset_idx}"
        # OJO: Aqui borramos toda la traza; incluida la del paso 1 (en lugar de ponerla a _ready). Porque posteriormente
        # vamos a generar la traza para este paso con estado '_in_progress'
        try:
            exec_query(db, query)
        except Exception as e:
            print(f"ERROR: can't delete traces : {str(e)}")

    def delete_results(config, sample, workflow):
        """ Remove existing results """

        # TODO: Implementar el borrado de los resultados downstream
        #       Actualmente solo borramos los resultados del step que vamos a lanzar. Deberiamos borrar tambien
        #       los resultados de los step posteriores
        ruta = root_dir(config, sample, workflow)
        shutil.rmtree(ruta, ignore_errors=True, onerror=None)

    def samples_to_be_processed(runid, sample, reset_idx):
        """
        Get the sample_id of samples to be processed. Samples in `_in_process` will not be processed

        :param runid: runid if set, None otherwise
        :param sample:  sample if set, None otherwise
        :return: tuple:
        """

        try:
            if not runid and not sample:
                # All samples
                if reset_idx is None:
                    step_idx = final_process_step()[1]
                    query = f"SELECT distinct(sample) FROM Tracks WHERE project = '{project_name}' " \
                            f"AND state IN('_ready', '_successful') AND current_step < {step_idx} " \
                            f"UNION SELECT distinct(sample) FROM Tracks WHERE  project = '{project_name}' " \
                            f"AND state IN('_aborted', '_failed')"
                else:
                    query = f"SELECT distinct(sample) FROM Tracker WHERE project = '{project_name}' " \
                            f"AND step_idx = {reset_idx} "

                return exec_query(db, query) if batch_size is None else exec_query(db, query)[:batch_size]

            elif not runid and sample:
                # Just this sample
                # query = "PENDIENTE DE IMPLELEMTACION!"
                # return exec_query(db, query)
                return [(sample,)]
            elif runid and not sample:
                # All samples of a Run
                query = "PENDIENTE DE IMPLELEMTACION!"
                return exec_query(db, query) if batch_size is None else exec_query(db, query)[:batch_size]
            else:
                # No sense to specify both runid and sample
                # TODO: manejar esta situacion en argparse
                print(f"ERROR: No sense to specify both {sample} and {runid}.")
                exit(1)
        except Exception as e:
            print(f"ERROR: can't read Tracks data : {str(e)}")
            exit(1)

    def next_workflow_to_run(sample_id, reset_idx):
        """ Get the next workflow to be run"""

        def get_info_from_run(sample_id):
            """ Given a sample get its top downstream execution state and step_idx"""
            query = f"""SELECT state, current_step FROM Tracks WHERE sample = '{sample_id}'"""
            try:
                return exec_query(db, query)[0]
            except Exception as e:
                print(f"ERROR: can't read Tracks info : {str(e)}")
                exit(1)

        if reset_idx:
            # Reset to `reset_idx` step
            idx = reset_idx
        else:
            state, step = get_info_from_run(sample_id)
            if state in ('_ready', '_aborted', '_failed', '_in_progress'):
                # Get the current step
                idx = step
            elif state == '_successful':
                # Get the next step
                idx = step + 1

        switcher = {
            1: ('consensus_reads', 1),
            2: ('variant_discovery', 2),
            3: ('somatic_calling', 3)
        }

        return switcher.get(idx, None)

    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')

    for sample in samples_to_be_processed(runid, sample, reset_idx):
        # TODO: Incluir logging

        sample_id = sample[0]
        workflow = next_workflow_to_run(sample_id, reset_idx)
        if reset_idx:
            # Delete info of previous executions
            delete_trace(config, sample_id, reset_idx)
            # delete_results(config, sample_id, workflow[0])
        if workflow:
            delete_results(config, sample_id, workflow[0])
            chip_bk.process(config, db, sample_id, workflow, flag_test)
        else:
            print(f"No next workflow to run for {sample_id}.")


def abort(config, runid, sample, flag_force=False):
    def samples():
        if sample is not None and runid is not None:
            # No sense to specify both runid and sample
            # TODO: manejar esta situacion en argparse
            print(f"ERROR: No sense to specify both {sample} and {runid}.")
            exit(1)
        elif sample is not None and runid is None:
            query = f"SELECT uuid  FROM Tracks WHERE state = '_in_progress' AND project = '{project_name}' " \
                    f"AND sample = '{sample}'"

        elif sample is None and runid is not None:
            query = f"SELECT uuid FROM Tracks  WHERE state = '_in_progress' " \
                    f"AND sample IN (SELECT sample FROM samples WHERE inclusion='enabled' AND project = '{project_name}' " \
                    f"AND runid = '{runid}')"
        else:
            query = f"""SELECT uuid  FROM Tracks WHERE state = '_in_progress' AND project = '{project_name}'"""

        try:
            result = exec_query(db, query)
        except Exception as e:
            print(f"ERROR: can't read sample data : {str(e)}")

        if len(result) == 0:
            print(f'ERROR: There is no samples \'_in_progress\'')
            exit(1)
        else:
            return result

    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')

    if sample is None and runid is None and not flag_force:
        print('ERROR: must use --force to apply all \'_in_progress\' samples')
        exit(1)
    else:
        for uuid in samples():
            chip_bk.abort(config, db, uuid[0])


def check(config, cleanup, flag_summary):
    def samples():
        query = f"SELECT sample, current_step, uuid FROM Tracks " \
                f"WHERE project = '{project_name}' and state = '_in_progress'"
        # TEST: return [('AA00293282', 1, '4a6979b3-d2b7-44a6-9228-c5341c1738bb')]
        try:
            return exec_query(db, query)
        except Exception as e:
            print(f"ERROR: can't read sample data : {str(e)}")
            exit(1)

    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')

    for sample, current_step, uuid in samples():
        chip_bk.check(config, db, project_name, sample, current_step, uuid)

    if flag_summary:
        pass
        # TODO: Implement 'summary'


def cohort(config, cohort, samples_file, flag_replace):
    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')

    if flag_replace:
        sql = f"DELETE FROM cohorts WHERE project = '{project_name}' AND cohort = '{cohort}'"
        db.cursor().execute(sql)

    try:
        with open(samples_file, 'r') as samples:
            for sample in samples:
                chip_bk.cohort(db, project_name, cohort, sample.strip())
        db.commit()
    except Exception as e:
        db.rollback()
        print(f"ERROR: {str(e)}")
    finally:
        db.close()


def report(config):
    chip_bk.report(config)


def candidates(config, cohort, flag_test=False):
    """ Generar listado de mutaciones CHIP"""

    def get_samples():
        """ Devuelve un array con el path al fichero de entrada correspondiente a cada muestra"""

        if cohort:
            query = f"SELECT DISTINCT(T.sample) FROM Tracks T INNER JOIN cohorts C ON T.project=C.project " \
                    f"AND T.sample=C.sample WHERE current_step = {last_wf_idx} " \
                    f"AND state ='_successful' AND project = '{project_name}' AND cohort = '{cohort}'"
        else:
            # Originalmente incluimos la muestras que estan en alguna cohorte; pero como esto no lo estan manejando aún
            # muy bien, lo que hago es obtener todas las muestra para las que se ha finalizado en procesamiento
            # query = "SELECT distinct(T.sample) FROM Tracks T INNER JOIN cohorts C ON T.sample=C.sample WHERE current_step = {} " \
            #        "AND state ='_successful'".format(final_process_step())
            query = f"SELECT distinct(T.sample) FROM Tracks T WHERE current_step = {last_wf_idx} " \
                    f"AND project = '{project_name}' AND state ='_successful'"
        try:
            return exec_query(db, query)
        except Exception as e:
            print(f"ERROR: can't read sample data : {str(e)}")
            exit(1)

    def get_path(sample):
        """
                Devuelve el path al fichero de la muestra. Comprueba que el fiochero existe; y si no existe
                devuelve un error.
            """

        path = f"{root_dir(config, sample, last_wf_name)}/_out/{sample}.vcf.gz"
        if not Path(path).exists():
            print(f'ERROR: Missing mutations for {sample} in cohort.')
            exit(1)
        else:
            return path

    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')

    last_wf_name, last_wf_idx = final_process_step()

    samples = [get_path(sample) for sample, in get_samples()]

    # samples = ["/data3/200831_JFuster_Somatic_PESA/ANALYSIS/AA00326716/somatic_calling/_out/AA00326716.vcf.gz",
    #           "/data3/200831_JFuster_Somatic_PESA/ANALYSIS/AA00326717/somatic_calling/_out/AA00326717.vcf.gz",
    #           "/data3/200831_JFuster_Somatic_PESA/ANALYSIS/AA00326718/somatic_calling/_out/AA00326718.vcf.gz"]
    if samples:
        cohort_name = cohort if cohort else 'all_{}'.format(project_name)
        chip_bk.candidates(config, db, cohort=cohort_name, samples=samples, flag_test=flag_test)
        db.commit()
    else:
        db.close()
        print(f"ERROR: there is no samples in cohort {cohort}.")
        exit(1)

    db.close()


def fingerprint(config, flag_test=False):
    """
        Genera un VCF con polimormismos para usarse como 'huella' para controlar
        la identidad de las muestras (por ejemplo, diferenmtes visitas del mismo individuos
    """

    def samples_to_be_processed():
        with open(samples_in, 'r') as samples_f:
            return samples_f.read().splitlines()

    # samples_in = '/home/jdelabarrera/DEV/check_RNA_SOM_PESA/Germline_calls/samples_to_call.tsv'
    samples_in = '/home/jdelabarrera/DEV/check_RNA_SOM_PESA/Germline_calls/second_round.tsv'

    for sample in samples_to_be_processed():

        if Path(f"/data3/200831_JFuster_Somatic_PESA/_consensus/{sample}.bam").exists():
            # TODO: Incluir logging
            chip_bk.fingerprint(config, sample, flag_test)
        else:
            print(f"Can't call germline variant call for {sample}.")


def dummy(config, bypass_vaf_filter, annotated_vcf, tsv_outfile):
    """

    :param config_file:
    :return:
    """

    chip_parsing_longitudinal.parse(config, annotated_vcf, tsv_outfile, bypass_vaf_filter)


def crossectional(config, bypass_vaf_filter, translation_filename, annotated_vcf, tsv_outfile):
    chip_crossectional.parse(config, annotated_vcf, tsv_outfile, translation_filename, bypass_vaf_filter)


def longitudinal(config, bypass_vaf_filter, translation_filename, annotated_vcf, prefix_outfile):
    chip_longitudinal.parse(config, annotated_vcf, prefix_outfile, translation_filename, bypass_vaf_filter)


def stat_cov(config, cohort, stat_outfile):
    # TODO: Checks if the cohort exists
    db = connect_db(config)
    chip_bk.stat_cov(config, db, cohort, stat_outfile)


def dummy_crossectional(config, bypass_vaf_filter, annotated_vcf, tsv_outfile):
    """

    :param config_file:
    :return:
    """

    crossectional.parse(config, annotated_vcf, tsv_outfile, bypass_vaf_filter)


def dummy_clean(config):
    """
        Copia cada CRAM a el directorio del proyecto correspondiente
    :param config:
    :return:
    """

    def samples_by_project(project):
        query = f"select distinct(sample) from CHIP_INSPECTOR.samples where project = '{project}' " \
                f"and inclusion = 'enabled'"
        return exec_query(db, query)

    db = connect_db(config)
    projects = ('HF_1', 'PESA_SOMATIC', 'Somatic_DP', 'Somatic_TB')
    orig = '/data3/200831_JFuster_Somatic_PESA/_consensus'

    for project in projects:
        dest_dir = f"{orig}/{project}"
        if not os.path.isdir(dest_dir):
            os.mkdir(dest_dir)
        for sample in samples_by_project(project):
            org_file = f"{orig}/{sample[0]}.cram"
            dest_file = f"{dest_dir}/{sample[0]}.cram"
            if Path(org_file).exists():
                shutil.move(org_file, dest_file)
    db.close()


def dummy_BAK(config, runid):
    """ Performs `dummy` operations. """

    def samples_to_be_processed(runid):
        # All samples of a Run
        query = f"""SELECT sample, current_step, uuid FROM 'Tracks' 
            WHERE state IN('_ready', '_aborted', '_failed', '_successful')
             AND sample IN (SELECT sample FROM samples WHERE inclusion='enabled' AND runid = '{runid}')"""
        return exec_query(query)

    for sample, current_stepm, uuid in samples_to_be_processed(runid):
        sample_id = sample[0]

        print(f"Borra {runid} - {sample}")
        # Borra resultados
        ruta = root_dir(config, sample, 'consensus_reads')
        shutil.rmtree(ruta, ignore_errors=True, onerror=None)
        # BORRAMOS la ejecucion de cromwell-executions
        dir = f"/data_PESA/cromwell-executions/consensus_reads/{uuid}"
        shutil.rmtree(dir, ignore_errors=True, onerror=None)

    # chip_bk.dummy(config)


def dummy_load_basestats(config):
    """ Load base statistics:
            - cobertura por gen (obviando chrX y chrY)
            - cobertura por muestra/gen"""

    def create_translation_table_from_file(translation_filename):
        res = dict()
        try:
            with open(translation_filename, 'r') as translation_file:
                for line in translation_file:
                    orig, dest = line.strip().split('\t')
                    res[orig] = dest
            return res
        except OSError as e:
            print(f"Unable to open {translation_filename}: {e}", file=sys.stderr)
            return dict()

    def samples():
        current_step = final_process_step()[1]
        query = f"SELECT sample FROM Tracks " \
                f"WHERE project = '{project_name}' and current_step = {current_step}"
        return exec_query(db, query)

    db = connect_db(config)
    project_name = config.get('ACCOUNTING', 'project_name')

    # Carga la tabla de conversion de nombre de muestra
    translation_filename = config.get('DIR', 'conversion_table')
    conversion_table = create_translation_table_from_file(translation_filename)

    for sample, in samples():
        chip_bk.dummy_load_basestats(config, db, project_name, sample, conversion_table)

    # print("sample\tAMELY\tBCOR\tBCORL1\tBEX2\tBEX4\tBRCC3\tDDX3X\tEIF1AY\tHSFX1\tKDM6A\tSMPX\tSRY\tZFY\tZRSR2\n")
    # for sample, in samples():
    #
    #     query = "select gene, cov_consensus from CHIP_INSPECTOR.sample_coverage_by_gene " \
    #         f"where sample = '{sample}' and gene in ('BCOR','BCORL1','BEX2','BEX4','BRCC3','DDX3X','HSFX1','KDM6A','SMPX','ZRSR2','AMELY'," \
    #         "'EIF1AY','SRY','ZFY') order by gene asc"
    #
    #     print(f"{sample}", end='')
    #     for gene, cov in exec_query(db, query):
    #         print(f"\t{cov}", end='')
    #     print(f"\n", end='')

    db.close()


def dummy_target_panel(config):
    """ Pruebas con la clase TargetPanel"""
    panel = tp.TargetPanel(config)

    genes = panel.genes()
    transcripts = panel.transcripts()

    k = 0


def dummy_check_somatic_dp(config):
    """ Para organizar los datos de todas las muestras de Somatic_DP voy a generar una tabla para identificarlas"""

    def get_sample_info():
        query = f"select sample , runid, experiment  from CHIP_INSPECTOR.samples s where project ='SOMATIC_DP'" \
                f"and sample = '{sample}'"
        return exec_query(db, query)

    db = connect_db(config)
    out = dict()
    base_dir = '/data3/20220825_DPascual_DeepSomatic_1/ANALYSIS'
    #sbase_dir = '/data3/200831_JFuster_Somatic_PESA/ANALYSIS'
    for root, directories, files in os.walk(base_dir, topdown=False, followlinks=True):
        items = root.split('/')
        if len(items) > 6:
            if root.endswith('variant_discovery/_out'):
                sample = items[4]
                if db_data := get_sample_info():
                    out[sample] = {
                        'runid': db_data[0][1],
                        'experiment': db_data[0][2],
                        'cram': 'N',
                        'gVCF': 'N',
                        'VCF': 'N'
                    }
                    for file in files:
                        if file.endswith('.cram'):
                            out[sample]['cram'] = 'Y'
                        elif file.endswith('.g.vcf.gz'):
                            out[sample]['gVCF'] = 'Y'
                        elif file.endswith('.vcf.gz'):
                            out[sample]['VCF'] = 'Y'
    print("sample\tRunID\tExperiment\tcram\tgVCF\tVCF")
    for sample, data in out.items():
        print(f"{sample}\t{data['runid']}\t{data['experiment']}\t{data['cram']}\t{data['gVCF']}\t{data['VCF']}")

