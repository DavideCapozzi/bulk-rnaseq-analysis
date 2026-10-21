# Configuration, explicit resource retrieval and capability resolution.


def write_if_changed(path, text):
    """Rewrite path only when its content differs: its mtime then re-triggers exactly the jobs that read it.

    Snakemake's params/code triggers depend on the metadata in .snakemake/, which is lost when the working
    folder changes; content-addressed inputs do not. Runs at parse time, so a dry run also refreshes them.
    """
    if os.path.exists(path):
        with open(path) as fh:
            if fh.read() == text:
                return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


# Resolved config (after schema defaults) and toolkit identity: a change in either re-runs the whole workflow,
# so the results always come from the code recorded in toolkit.json.
write_if_changed(RESOLVED, yaml.safe_dump({**config, "_config_dir": CONFIG_DIR}, sort_keys=False))
write_if_changed(TOOLKIT_JSON, json.dumps(TOOLKIT, indent=2) + "\n")


rule fetch_resources:
    # Not an input of any analysis rule: run it explicitly. Analysis steps only read the cache and verify md5.
    input:
        config=RESOLVED,
        toolkit=TOOLKIT_JSON,
    output:
        os.path.join(CACHE, "resources_manifest.tsv"),
    log:
        os.path.join(PROV, "logs", "fetch_resources.log"),
    shell:
        "{RSCRIPT:q} {STEP:q} fetch_resources {input.config} {REPO} > {log} 2>&1"


checkpoint capabilities:
    input:
        config=RESOLVED,
        toolkit=TOOLKIT_JSON,
    output:
        os.path.join(PROV, "capabilities.json"),
    log:
        os.path.join(PROV, "logs", "capabilities.log"),
    shell:
        "{RSCRIPT:q} {STEP:q} capabilities {input.config} {REPO} > {log} 2>&1"


def load_capabilities():
    with open(checkpoints.capabilities.get().output[0]) as fh:
        return json.load(fh)


def done_file(step, caps):
    if step == "ingest":
        return os.path.join(OUT, "01_ingest", "step.done")
    row = next(r for r in caps if r["module"] == step)
    return os.path.join(OUT, row["dir"], "step.done")


rule plan:
    # Target of `bulkrnaseq plan`: prints which modules run and why the others are skipped. No output, so it
    # always runs; nothing else depends on it.
    input:
        os.path.join(PROV, "capabilities.json"),
    run:
        with open(input[0]) as fh:
            caps = json.load(fh)
        width = max(len(r["dir"]) for r in caps)
        print(f"\nPlan for {config['project']['name']} (toolkit {TOOLKIT['version']}, dev mode: {TOOLKIT['dev_mode']})")
        for r in caps:
            print(f"  {r['dir']:<{width}}  {r['status']:<8}  {r['reason'] or ''}")
