# Nextflow pipeline conventions at Nexomis

## Workflow file structure
```
.
└── <worflow_name>/
    ├── conf/
    │   ├-- ext.conf
    │   ├-- publish.conf
    ├── modules/
    │   ├-- config/
    │   │   ├── process/
    │   │   │   └── labels.config
    │   │   └── profiles.config
    │   ├-- process/
    │   │   ├── <process_name_1>/
    │   │   │   └── main.tf
    │   │   └── <process_name_2>/
    │   │       └── main.tf
    │   └-- subworflow/
    │       └── <subworflow_name_1>/
    │           └── main.tf
    ├── main.nf
    ├── nextflow.config
    ├── nextflow_schema.json
    └── Readme.md
```

## Help
The `help` will be managed by `plugin/nf-schema`:
  - Include in `main.nf`
  - Create the `nextflow_schema.json` file according to the parameter values in the `params{}` section in `nextflow.config`

## Config
`nextflow.config` will include the following config files:
  - `modules/config/process/labels.config` for managing resources assigned to processes/labels
  - `modules/config/profiles.config`
  - `conf/publish.config` to manage the publication of processes and subworkflows
**NOTE:** configuration parameters can be overridden by the `-c` option in Nextflow.

## Modules

Except in very rare cases where no existing process or subworkflow calls are made, and where no creation of processes or subworkflows is relevant, the following conventions should be followed.

### Structure of `modules/
Since `include` works only at the repository level, modules will be linked in the `modules` folder (to avoid file duplication and prevent incompatibility issues during module updates).

```sh
mkdir modules
git submodule add https://github.com/nexomis/nf-subworkflows.git modules/subworkflows
git submodule add https://github.com/nexomis/nf-config.git modules/config
git submodule add https://github.com/nexomis/nf-process.git modules/process
```
**NOTE:** submodules are linked to a repository by their hash reference. For them to be updated, they need to be pulled/pushed again.

For each modularized element (subworkflow or subprocess), place it in a folder named after the element containing a `main.nf`, which will be sourced for calls.

### Continue modularity: reusable, and easily maintainable

  - Tag processes with labels present in `modules/config/process/labels.config` and update if necessary.

  - To easily adapt (without multiplying module versions) the publication of `processes` and `subworkflows` (not the case for `workflows`) specifically to each workflow, centralize publication operations in the `conf/publish.conf` file.  
**Note:** to target a specific call of a process or subworkflow they must be imported with an unique name. In all cases a process (or subworkflow) can't be call twices with the same name.

  - Do not directly call parameters (`params`) in processes and sub-workflows: processes and sub-workflows should not directly depend on global parameters. Use channels to pass parameter values in workflows: aggregate parameters into channels at the main workflow level, then pass these channels to processes and sub-workflows as inputs.

  - In order to pass arguments to specific functions use `conf/ext.conf`, those arguments can then be access directly in the subworkflows or processes by calling `task.ext.args`.

## Notations

  - **Reads:** in every workflow, reads will be defined as a tuple with the following structure: `(sample_name, [path(s)_to_file(s)])`

  - **Typography:** process and workflow names will be written in capital letters: `workflow MY_WF {}`, `process MY_PROCESS {}`.


