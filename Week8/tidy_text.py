import sys, textwrap, re

def tidy_text(file_name):

    f = open(file_name, 'r')

    if f is None:
        sys.stderr("Can't open text file" , file_name, " for reading")

    lines = f.readlines()
    f.close()

    contents = " ".join(lines)

    # First, convert all uppercase text to lowercase
    contents = contents.lower()

    # Second, convert dashes (--) into spaces
    contents = re.sub("--", " ", contents)

    # ...and then join hyphenated words into single words
    contents = re.sub(r"-\s*", "", contents)

    # Third, convert all whitespace characters into spaces
    contents = re.sub(r"\s", " ", contents)

    # ...and then remove anything that is not a letter or space
    contents = re.sub(r"[^a-z ]", "", contents)

    # Fourth, convert any run of multiple spaces into a single space
    contents = re.sub(r"\s+", " ", contents)

    # Finally, break lines appropriately with newlines; if you would
    # prefer the output to be on one very long line (with no newline
    # characters; only spaces) because you think this would make your
    # processing easier, just comment out these two lines

    content_data = textwrap.wrap(contents, 72)

    filename = "tidy." + file_name
    f = open(filename, "w")

    if f is None:
        sys.stderr("output file not found: tidy." + file_name)
    else:
        # ...and output the results
        for line in content_data:
            f.write(line + " \n")
        f.close()
        print("Wrote %s" % filename)


if __name__ == '__main__':
    # change the parameters if necessary
    tidy_text('paradise.lost.txt')
    tidy_text('heart.of.darkness.txt')
