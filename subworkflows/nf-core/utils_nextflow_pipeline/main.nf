//
// Subworkflow with functionality that may be useful for any Nextflow pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTILS_NEXTFLOW_PIPELINE {
    take:
    print_version        // boolean: print version
    dump_parameters      // boolean: dump parameters
    outdir               //    path: base directory used to publish pipeline results
    check_conda_channels // boolean: check conda channels
    params_map           //     map: copy of the params

    main:

    //
    // Print workflow version and exit on --version
    //
    if (print_version) {
        log.info("${workflow.manifest.name} ${getWorkflowVersion()}")
        System.exit(0)
    }

    //
    // Dump pipeline parameters to a JSON file
    //
    if (dump_parameters && outdir) {
        dumpParametersToJSON(outdir, params_map)
    }

    //
    // When running with Conda, warn if channels have not been set-up appropriately
    //
    if (check_conda_channels) {
        checkCondaChannels()
    }

    emit:
    dummy_emit = true
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get commit hash
//
def getGitHash() {
    if (workflow.commitId) {
        return workflow.commitId
    }
    def sout = new StringBuffer()
    def serr = new StringBuffer()
    def proc = "git rev-parse HEAD".execute(null, file(workflow.projectDir).toFile())
    proc.waitForProcessOutput(sout, serr)
    return proc.exitValue() ? null : "${sout}".trim()
}

//
// Get current branch
//
def getGitRevision() {
    if (workflow.revision) {
        return workflow.revision
    }
    def sout = new StringBuffer()
    def serr = new StringBuffer()
    def proc = "git rev-parse --abbrev-ref HEAD".execute(null, file(workflow.projectDir).toFile())
    proc.waitForProcessOutput(sout, serr)
    return proc.exitValue() ? null : "${sout}".trim()
}

//
// Generate version string
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    def commitId = getGitHash()
    if (commitId) {
        def git_shortsha = commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    def revision = getGitRevision()
    if (revision) {
        version_string += "@${revision}"
    }

    return version_string
}

//
// Dump pipeline parameters to a JSON file
//
def dumpParametersToJSON(outdir, params_map) {
    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    def filename = "params_${timestamp}.json"
    def temp_pf = new File(workflow.launchDir.toString(), ".${filename}")
    def jsonStr = groovy.json.JsonOutput.toJson(params_map)
    temp_pf.text = groovy.json.JsonOutput.prettyPrint(jsonStr)

    nextflow.extension.FilesEx.copyTo(temp_pf.toPath(), "${outdir}/pipeline_info/params_${timestamp}.json")
    temp_pf.delete()
}

//
// When running with -profile conda, warn if channels have not been set-up appropriately
//
def checkCondaChannels() {
    def parser = new org.yaml.snakeyaml.Yaml()
    def channels = []
    try {
        def config = parser.load("conda config --show channels".execute().text)
        channels = config.channels
    }
    catch (NullPointerException _e) {
        log.warn("Could not verify conda channel configuration.")
        return null
    }
    catch (IOException _e) {
        log.warn("Could not verify conda channel configuration.")
        return null
    }

    // Check that all channels are present
    // This channel list is ordered by required channel priority.
    def required_channels_in_order = ['conda-forge', 'bioconda']
    def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

    // Check that they are in the right order
    def channel_priority_violation = required_channels_in_order
        .withIndex()
        .any { channel, index ->
            index < required_channels_in_order.size() - 1 && !(channels.indexOf(channel) < channels.indexOf(required_channels_in_order[index + 1]))
        }

    if (channels_missing | channel_priority_violation) {
        log.warn(
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  There is a problem with your Conda configuration!\n\n" + "  You will need to set-up the conda-forge and bioconda channels correctly.\n" + "  Please refer to https://bioconda.github.io/\n" + "  The observed channel order is \n" + "  ${channels}\n" + "  but the following channel order is required:\n" + "  ${required_channels_in_order}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        )
    }
}
