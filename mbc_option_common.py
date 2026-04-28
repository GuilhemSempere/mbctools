def prepare_with_previous_params(core):
    core.prevent()
    core.prev_param(None)


def show_missing_loci_and_return(core, message, return_menu_callback):
    print(core.errorStyle + message + core.normalStyle)
    input("\nPress ENTER to continue ")
    return_menu_callback(core)
