from parsing.function_patcher import patch_function

# snakemake function execution
paths = {
    'bhsadata': snakemake.input.bhsadata,
    'editfuncts': snakemake.input.editfuncts,
    'patched': snakemake.output.patched,
}

patch_function(paths)
