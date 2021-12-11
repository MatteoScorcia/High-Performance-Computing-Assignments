import os


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


if __name__ == '__main__':
    directory = os.path.realpath("assigment1/section2/csv")
    for filename in os.listdir(directory):
        full_path = os.path.join(directory, filename)

        with open(full_path, "r") as f:
            contents = f.readlines()

        contents.insert(6, ",")

        with open(full_path, "w") as f:
            contents = "".join(contents)
            f.write(contents)
