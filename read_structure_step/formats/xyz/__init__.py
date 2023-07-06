from read_structure_step.formats.registries import register_format_checker
from . import xyz  # noqa: F401


@register_format_checker(".xyz")
def check_format(file_name):
    with open(file_name, "r") as f:
        lines = f.read().splitlines()

        if len(lines) < 3:
            return False

        # Standard XYZ
        len1 = len(lines[0].split())
        fields3 = lines[2].split()
        len3 = len(fields3)
        if len1 == 1 and len3 == 4:
            try:
                int(lines[0].split())
            except Exception:
                pass
            else:
                return True

        # Minnesota variant
        if len(lines) > 3 and len(fields3) == 2:
            try:
                int(fields3[0])
                int(fields3[1])
            except Exception:
                pass
            else:
                if len(lines[3].split()):
                    return True

        return False
