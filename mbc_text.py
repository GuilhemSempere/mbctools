def dos2unix(file_path):
    with open(file_path, "rb") as file:
        content = file.read().decode("utf-8").replace("\r\n", "\n")
    with open(file_path, "wb") as file:
        file.write(content.encode("utf-8"))


def replaceInFile(filename, replacements):
    if len(replacements) == 0:
        return
    with open(filename, "r") as file:
        content = file.read()

    for old_text, new_text in replacements:
        content = content.replace(old_text, new_text)

    with open(filename, "w") as file:
        file.write(content)
