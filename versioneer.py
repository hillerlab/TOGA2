#!/usr/bin/env python3

"""Version update routine"""

import click
import in_place
import sys

CHANGELOG: str = "changelog.md"
README: str = "README.md"
VERSION_FILE: str = "__version__.py"
VERSION_TEMPLATE: str = "__version__ = \"{}\""
CHANGELOG_LINK_LINE: str = (
    "\nFor the full list of code changes, see "
    "[changelog.md](https://github.com/hillerlab/TOGA2/blob/main/changelog.md) .\n"
)

def remove_version_suffix(version: str) -> str:
    """
    Removes the letter suffix from the version name

    Args:
        version: version identifier
    
    Returns:
        Version identifier without letter suffices
    """
    return "".join((x for x in version if not x.isalpha()))

@click.command()
@click.argument("version", type=str, metavar="VERSION_NAME")
def versioneer(version: str) -> None:
    """Updates version and README.md"""
    version = version.strip()
    if not version[0].isdigit():
        raise ValueError("Version identifier must start with a digit")
    version = "v" + version
    version_main: str = remove_version_suffix(version)

    ## update __version__.py
    with open(VERSION_FILE, "w") as h:
        h.write(VERSION_TEMPLATE.format(version))

    ## fetch the recent update description from changelog.md
    changelog_lines: str = ""
    with open(CHANGELOG, "r") as h:
        for line in h:
            _line: str = line.strip()
            if line.startswith("#") and _line.endswith(version) or _line.endswith(version_main):
                changelog_lines += line
                continue
            if line.startswith("#") and changelog_lines or line.startswith("* Minor changes"):
                break
            if changelog_lines:
                changelog_lines += line
    if not changelog_lines:
        click.echo("WARNING: No changelog update found")
        sys.exit(0)
    changelog_lines += CHANGELOG_LINK_LINE

    ## replace the previous change description in README.md with a recent 
    anchor_found: bool = False
    with in_place.InPlace(README) as h:
        for line in h:
            if line.startswith("## Changelog"):
                h.write(line)
                anchor_found = True
                continue
            if anchor_found:
                if line.startswith("## "):
                    h.write("#" + changelog_lines + '\n')
                    anchor_found = False
                    h.write(line)
                else:
                    continue
            else:
                h.write(line)
            #     if not line.startswith("#"):
            #         continue
            #     h.write("#" + changelog_lines + '\n')
            #     anchor_found = False
            # else:
            #     h.write(line)

if __name__ == '__main__':
    versioneer()