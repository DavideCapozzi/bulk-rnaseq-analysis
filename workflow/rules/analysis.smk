# Ingest, analysis modules (one generic rule driven by the registry), provenance and report.


rule ingest:
    input:
        config=RESOLVED,
        capabilities=os.path.join(PROV, "capabilities.json"),
    output:
        os.path.join(OUT, "01_ingest", "step.done"),
    log:
        os.path.join(OUT, "01_ingest", "logs", "ingest.log"),
    shell:
        "{RSCRIPT:q} {STEP:q} ingest {input.config} {REPO} > {log} 2>&1"


def module_inputs(wildcards):
    caps = load_capabilities()
    row = next(r for r in caps if r["dir"] == wildcards.moddir)
    deps = [d for d in row["depends_on"].split(",") if d]
    return [RESOLVED, TOOLKIT_JSON] + [done_file(d, caps) for d in deps]


rule module:
    input:
        module_inputs,
    output:
        os.path.join(OUT, "{moddir}", "step.done"),
    log:
        os.path.join(OUT, "{moddir}", "logs", "module.log"),
    wildcard_constraints:
        moddir=r"(0[2-9]|1[0-9])_[a-z_]+",
    params:
        module=lambda wildcards: wildcards.moddir[3:],
    shell:
        "{RSCRIPT:q} {STEP:q} module {RESOLVED} {REPO} {params.module} > {log} 2>&1"


def run_modules_done(wildcards):
    caps = load_capabilities()
    return [os.path.join(OUT, "01_ingest", "step.done")] + [
        os.path.join(OUT, r["dir"], "step.done") for r in caps if r["status"] == "run"
    ]


rule provenance:
    input:
        run_modules_done,
    output:
        os.path.join(PROV, "run_manifest.json"),
    log:
        os.path.join(PROV, "logs", "provenance.log"),
    shell:
        "{RSCRIPT:q} {STEP:q} provenance {RESOLVED} {REPO} > {log} 2>&1"


rule report:
    input:
        os.path.join(PROV, "run_manifest.json"),
    output:
        os.path.join(OUT, "report", "summary.md"),
    log:
        os.path.join(PROV, "logs", "report.log"),
    shell:
        "{RSCRIPT:q} {STEP:q} report {RESOLVED} {REPO} > {log} 2>&1"
