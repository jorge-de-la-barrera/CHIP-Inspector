#!/usr/bin/env python3

import argparse
import chip_api
import configparser

# TODO: usar logging para manejar la salida al usuario


main_parser = argparse.ArgumentParser()

main_parser.add_argument('--config',
                         dest='config_file',
                         required=False,
                         default='../config.ini',
                         type=str.strip,
                         help='Path to the config file')
main_parser.add_argument('-P',
                    '--project',
                    dest='project',
                    required=False,
                    type=str.strip,
                    help='Project (by default that is specified in the .ini config file)')

subparsers = main_parser.add_subparsers(
    help='sub-command help', dest='command')

# commands
init = subparsers.add_parser(
    'init', help='init help', description='Initialize project.')
load = subparsers.add_parser(
    'load', help='load help', description='Load uBAM sequence files to the project.')
drop = subparsers.add_parser(
    'drop', help='drop help', description='Discards samples from the project.')
enable = subparsers.add_parser(
    'enable', help='enable help', description='Set sample available for processing.')
disable = subparsers.add_parser(
    'disable', help='disable help', description='Set sample unavailable for processing.')
createPanel = subparsers.add_parser(
    'createPanel', help='createPanel help', description='Create a new target panel.')
removePanel = subparsers.add_parser(
    'removePanel', help='removePanel help', description='Remove a target panel.')
process = subparsers.add_parser(
    'process', help='process help', description='Process samples.')
abort = subparsers.add_parser(
    'abort', help='abort help', description='Request Cromwell to abort a running workflow.')
check = subparsers.add_parser(
    'check', help='check help', description='Checks and updates status of processing.')
report = subparsers.add_parser(
    'report', help='report help', description='Shows a reports with results.')
candidates = subparsers.add_parser(
    'candidates', help='candidates help', description='Generation of CHIP mutations candidates.')
cohort = subparsers.add_parser(
    'cohort', help='cohort help', description='Manage cohorts within a project.')
fingerprint = subparsers.add_parser(
    'fingerprint', help='fingerprint help', description='Create a germinal call for QC.')
crossectional = subparsers.add_parser(
    'crossectional', help='crossectional', description='Create list of candidate mutation for a crossectional study.')
longitudinal = subparsers.add_parser(
    'longitudinal', help='longitudinal', description='Create list of candidate mutation (candidates) and passenger '
                                                     '(synonymous) for a longitudinal study.')
stat_cov = subparsers.add_parser(
    'statCov', help='statCov', description='Calculate the coverage statistics for a cohort.')
dummy = subparsers.add_parser(
    'dummy', help='dummy help', description='Performs `dummy` operations.')

# 'init' arguments
init.add_argument('--force',
                  dest='flag_force',
                  action='store_true',
                  required=True,
                  default=False,
                  help='To avoid rewrite existing data by error the --force option must be used.')
# 'load' arguments
load.add_argument('--samplesheet',
                  dest='samplesheet',
                  type=str.strip,
                  required=True,
                  help='Samplesheet of the Run to be loaded.')
load.add_argument('--project-to-import',
                  dest='samplesheet_project',
                  type=str.strip,
                  required=True,
                  help='Project of the samplesheet to be loaded.')
load.add_argument('-R',
                  '--runid',
                  dest='runid',
                  type=str.strip,
                  required=True,
                  help='ID of the Run to be loaded.')
load.add_argument('--skip-missing-samples',
                  dest='flag_skip_missing_samples',
                  action='store_true',
                  required=False,
                  default=False,
                  help='Continue processing samples if the sequences of a sample is missing.')


# 'drop' arguments
drop.add_argument('-R',
                  '--runid',
                  dest='runid',
                  required=False,
                  type=str.strip,
                  default=None,
                  help='Discard samples of the given <runid>.')
drop.add_argument('-S',
                  '--sample',
                  dest='sample',
                  required=False,
                  type=str.strip,
                  default=None,
                  help='Discard sample <sample>')
# 'enable' arguments
enable.add_argument('-R',
                    '--runid',
                    dest='runid',
                    required=False,
                    type=str.strip,
                    default=None,
                    help='Enable samples of the given <runid>.')
enable.add_argument('-S',
                    '--sample',
                    dest='sample',
                    required=False,
                    type=str.strip,
                    default=None,
                    help='Enable sample <sample>')
# 'disable' arguments
disable.add_argument('-R',
                     '--runid',
                     dest='runid',
                     required=False,
                     type=str.strip,
                     default=None,
                     help='Disable samples of the given <runid>.')
disable.add_argument('-S',
                     '--sample',
                     dest='sample',
                     required=False,
                     type=str.strip,
                     default=None,
                     help='Disable sample <sample>')
# createPanel arguments
createPanel.add_argument('-N',
                         '--panel_name',
                         dest='panel_name',
                         required=True,
                         type=str.strip,
                         default=None,
                         help='Name of the panel to be created.')
createPanel.add_argument('-F',
                         '--file',
                         dest='panel_file',
                         required=False,
                         type=str.strip,
                         default=None,
                         help='Path to the file with the panel definition (columns: 1-Gene symbol, 2-RefSeq feature).')
#createPanel.add_argument('--force',
#                         dest='flag_force',
#                         required=False,
#                         type='store_true',
#                         default=False,
#                        help='To avoid rewrite existing data by error the --force option must be used.')
#removePanel arguments
removePanel.add_argument('-N',
                         '--panel_name',
                         required=True,
                         type=str.strip,
                         default=None,
                         help='Name of the panel to be removed.')
# 'process' arguments
process.add_argument('-R',
                     '--runid',
                     dest='runid',
                     required=False,
                     type=str.strip,
                     default=None,
                     help='Samples of the given <runid> run to be processed.')
process.add_argument('-S',
                     '--sample',
                     dest='sample',
                     required=False,
                     type=str.strip,
                     default=None,
                     help='Sample <sample> to be processed')
process.add_argument('--reset',
                     dest='reset_idx',
                     required=False,
                     type=int,
                     default=None,
                     help='Re-processed from the workflow specified (1 from the begining). Only affects to samples'
                          'not `_in_progress`.')
process.add_argument('--test',
                     dest='flag_test',
                     required=False,
                     action='store_true',
                     default=False,
                     help='Submits workflows to the \'test\' queue.')
process.add_argument('--batch_size',
                     dest='batch_size',
                     required=False,
                     type=int,
                     default=None,
                     help='Number of samples to be processed in a batch.')

# 'abort' arguments
abort.add_argument('-R',
                   '--runid',
                   dest='runid',
                   required=False,
                   type=str.strip,
                   default=None,
                   help='Samples of the given <runid> run to be processed.')
abort.add_argument('-S',
                   '--sample',
                   dest='sample',
                   required=False,
                   type=str.strip,
                   default=None,
                   help='Sample <sample> to be processed')
abort.add_argument('--force',
                   dest='flag_force',
                   action='store_true',
                   required=False,
                   default=False,
                   help='To avoid rewrite existing data by error the --force option must be used.')
# 'check' argument
check.add_argument('--summary',
                   dest='flag_summary',
                   action='store_true',
                   default=False,
                   help='Shows detailed information about the processing status.')
check.add_argument('--cleanup',
                   dest='cleanup',
                   action='store_true',
                   default=False,
                   help='Cleans transient data.')
# 'cohort' arguments
cohort.add_argument('-C',
                    '--cohort',
                    dest='cohort',
                    required=True,
                    type=str.strip,
                    default=None,
                    help='Cohort for which samples will be assigned.')
cohort.add_argument('--samples',
                    dest='samples_file',
                    type=str.strip,
                    required=True,
                    help='List of samples.')
cohort.add_argument('--replace',
                    dest='flag_replace',
                    action='store_true',
                    default=False,
                    help='Delete previous data.')
# 'candidates' arguments
candidates.add_argument('-C',
                        '--cohort',
                        dest='cohort',
                        required=False,
                        type=str.strip,
                        default=None,
                        help='Cohort for which CHIP mutations will be identified.')
candidates.add_argument('--test',
                        dest='flag_test',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Submits workflows to the \'test\' queue.')
# 'crossectional' arguments
crossectional.add_argument('--bypass_vaf_filter',
                   dest='bypass_vaf_filter',
                   required=False,
                   action='store_true',
                   default=False,
                   help='Bypass minimum VAF filter if set to TRUE.')
crossectional.add_argument('--translation-sample-names',
                   dest='translation_filename',
                   type=str.strip,
                   required=False,
                   default=None,
                   help='File for sample names translation (col 1: old name, col 2: new name).')
crossectional.add_argument('-V', '--vcf',
                   dest='annotated_vcf',
                   type=str.strip,
                   required=True,
                   help='VCF multisample of the cohort.')
crossectional.add_argument('-O', '--output',
                   dest='tsv_outfile',
                   type=str.strip,
                   required=True,
                   help='File with the ouput table.')
# 'longitudinal' arguments
longitudinal.add_argument('--bypass_vaf_filter',
                   dest='bypass_vaf_filter',
                   required=False,
                   action='store_true',
                   default=False,
                   help='Bypass minimum VAF filter if set to TRUE.')
longitudinal.add_argument('--translation-sample-names',
                   dest='translation_filename',
                   type=str.strip,
                   required=False,
                   default=None,
                   help='File for sample names translation (col 1: old name, col 2: new name).')
longitudinal.add_argument('-O', '--prefix_outfile',
                   dest='prefix_outfile',
                   type=str.strip,
                   required=False,
                   default=None,
                   help='Prefix for the out files (CHIP candidates and synonymous mutations')
longitudinal.add_argument('-V', '--vcf',
                   dest='annotated_vcf',
                   type=str.strip,
                   required=True,
                   help='VCF multisample of the cohort.')

stat_cov.add_argument('-C', '--cohort',
                   dest='cohort',
                   type=str.strip,
                   required=True,
                   help='Cohort for which coverage statistics will be computed.')
stat_cov.add_argument('-O', '--output',
                   dest='stat_outfile',
                   type=str.strip,
                   required=True,
                   help='File with the ouput coverage statistics (xlsx).')
# Dummies

# 'dummy' arguments
# dummy.add_argument('--bypass_vaf_filter',
#                    dest='bypass_vaf_filter',
#                    required=False,
#                    action='store_true',
#                    default=False,
#                    help='Bypass minimum VAF filter if set to TRUE.')
#
# dummy.add_argument('-V', '--vcf',
#                    dest='annotated_vcf',
#                    type=str.strip,
#                    required=True,
#                    help='VCF multisample of the cohort.')
#
# dummy.add_argument('-O', '--output',
#                    dest='tsv_outfile',
#                    type=str.strip,
#                    required=True,
#                    help='File with the ouput table.')

if __name__ == '__main__':

    commands = {
        "init": chip_api.init,
        "load": chip_api.load,
        "drop": chip_api.drop,
        "enable": chip_api.enable,
        "disable": chip_api.disable,
        "process": chip_api.process,
        "abort": chip_api.abort,
        "check": chip_api.check,
        "cohort": chip_api.cohort,
        "report": chip_api.report,
        "candidates": chip_api.candidates,
        "fingerprint": chip_api.fingerprint,
        "crossectional": chip_api.crossectional,
        "longitudinal": chip_api.longitudinal,
        "createPanel": chip_api.create_panel,
        "removePanel": chip_api.remove_panel,
        "statCov": chip_api.stat_cov,
        #"dummy": chip_api.dummy_crossectional,
        #"dummy": chip_api.dummy_target_panel
        "dummy": chip_api.dummy_check_somatic_dp
    }


    def error():
        print(f"Command <{args['command']}> doesn't exist.")


    args = vars(main_parser.parse_args())

    # Call the command
    func = commands.get(args['command'], error)
    del args['command']

    try:

        # Lee la configuracion y sobreescribe los parametros que se han pasado por parametro
        config = configparser.ConfigParser(interpolation=None)
        config.read(args['config_file'])
        args['config'] = config
        del args['config_file']

        if args.get('project', False):
            config.set('ACCOUNTING', 'project_name', args.get('project'))
        del args['project']

        func(**args)
    except BaseException as error:
        print(error)
