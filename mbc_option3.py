from mbc_option_common import prepare_with_previous_params
from mbc_sequences import concat_sequences_by_locus


def concat_3(core, loci=None):
    """Compiles sample sequences by locus in a unique file (option 3)."""
    return concat_sequences_by_locus(
        loci=loci,
        lociPEs=core.lociPEs,
        lociSEs=core.lociSEs,
        samples=core.samples,
        current_dir=core.current_dir,
        fileSep=core.fileSep,
        promptUser=core.promptUser,
        main=core.main,
        derep_several_fasta_files=core.derepSeveralFastaFiles,
        warningStyle=core.warningStyle,
        errorStyle=core.errorStyle,
        successStyle=core.successStyle,
        normalStyle=core.normalStyle,
        promptStyle=core.promptStyle,
    )


def main_menu3(core):
    """Displays submenu 3."""
    prepare_with_previous_params(core)
    core.os.system("cls" if core.winOS else "clear")
    print(
        core.titleStyle
        + "\n--- MENU 3: GENERATION OF A UNIQUE SEQUENCE FILE FOR EACH LOCUS (comprising all samples' data) ---"
        + core.normalStyle
        + "\n\n3a -> Process all loci at once\n\n"
        "3b -> Process loci one by one\n"
    )

    core.rmenu = core.promptUser(
        "Please select an option among those listed above",
        None,
        ["3a", "3b", "back", "home", "exit"],
        1,
        core.main,
        "",
    )
    concat_3(core, core.lociPEs + list(set(core.lociSEs) - set(core.lociPEs)) if core.rmenu == "3a" else None)
    main_menu3(core)


def menu3(core):
    """Runs option 3."""
    prepare_with_previous_params(core)
    concat_3(core)
    core.rerun(None)
