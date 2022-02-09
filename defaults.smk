## Load default config values
# - loading type errors like 'list indices must be integers or slices, not str'
#   probably means the configfile given has a list where defaults were a dict.
#   eg. config['chromosomes'] as list vs config['chromosomes']['hg38'] as list

# snakefmt cannot understand 'configfile:' inside an if statement and will fail on this script
if "chromosome_defaults" not in config:
    for ref in ["hg19", "hg38"]:
        def_ref_yaml = os.path.join(workflow.basedir, f"resources/config/genome_ref.{ref}.yaml")
        logger.warning(f"Loading reference defaults from: {def_ref_yaml}")
        configfile: def_ref_yaml
if "bams" not in config:
    logger.error(f"No patient bams given.  Using demo test case.")
    configfile: os.path.join(workflow.basedir, "resources/config/default_case.yaml")

configfile: os.path.join(workflow.basedir, "resources/config/default_run_params.yaml")
