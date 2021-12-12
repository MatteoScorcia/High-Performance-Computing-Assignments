import os


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


if __name__ == '__main__':
    directory = "1csv"
    for filename in os.listdir(directory):
        full_path = os.path.join(directory, filename)

        with open(full_path, "r") as f:
            lines = f.readlines()
            lines[0] = lines[0].replace(' ', ',')

        with open(full_path, "w") as f:
            f.writelines(lines)
