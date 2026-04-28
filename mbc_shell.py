import platform

winOS = "Windows" in platform.uname()


def start_log_redirect(filepath):
    """Starts redirecting log messages."""
    if not winOS:
        return "exec 3>&1 4>&2 >" + filepath + " 2>&1\n"
    return 'Clear-Content -Path "' + filepath + '" -Force -ErrorAction SilentlyContinue\n&{\n'


def end_log_redirect(filepath):
    """Stops redirecting log messages."""
    if not winOS:
        return "exec 1>&3 2>&4\n"
    return '} 2>&1 | ForEach-Object { if ($_ -is [System.Management.Automation.ErrorRecord]) { $_.Exception.Message } else { $_ } } | Out-File -FilePath "' + filepath + '" -Append\n'


def main_stream_message(message):
    """Displays a main stream message on the console."""
    if not winOS:
        return "printf \"" + message + "\" >&3\n"
    return "Write-Host -NoNewline \"" + message + "\"\n"


def logFileMessage(msg):
    """Returns shell code to append a section marker in log files."""
    if not winOS:
        return f"echo\necho\necho {msg}\necho " + ("-" * len(msg)) + "\necho\n"
    return f'Write-Output ""\nWrite-Output ""\nWrite-Output "{msg}"\nWrite-Output "' + ("-" * len(msg)) + '"\nWrite-Output ""\n'
