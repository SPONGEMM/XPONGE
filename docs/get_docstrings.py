import os
from pathlib import Path
import re

class Function:
    def __init__(self, name):
        self.name = name
        self.docstring = ""

class Class:
    def __init__(self, name):
        self.name = name
        self.docstring = ""
        self.functions = []

class Module:
    modules = {}
    def __init__(self, name, path):
        self.modules[name] = self
        self.path = path
        self.name = name
        self.docstring = ""
        self.classes = []
        self.functions = []
        self.submodules = []
        self.subpackages = []

def transform_parm_and_return(input_string):
    lines = input_string.split('\n')

    usual_lines = []
    return_lines = []

    parameters_lines = []


    for line in lines:
        match1 = re.match(r':param\s+([^:]+):\s*(.*)', line)
        match2 = re.match(r':return:\s*(.*)', line)
        if match1:
            if not parameters_lines:
                parameters_lines.append("\n\nparameters | explanation")
                parameters_lines.append("--- | ---")
            param_name = match1.group(1).strip()
            param_description = match1.group(2).strip()
            parameters_lines.append(f"{param_name} | {param_description}")
        elif match2:
            return_name = match2.group(1).strip()
            if not return_name:
                return_name = None
            return_lines.append(f"""\n\nreturn | {return_name}\n---|---""")
        else:
            usual_lines.append(line)

    result_string = '\n'.join([lines[0]] + parameters_lines + return_lines)

    return result_string

def process_docstring(string):
    pattern = re.compile(r'\\\n', re.DOTALL)
    string = re.sub(pattern, '', string)
    pattern = re.compile(r'\.\. TIP::\s*\n(.*?)\s*\n\s*\n', re.DOTALL)
    string = re.sub(pattern, lambda match: '> ' + match.group(1).replace('\n', '\n> ') + '\n{.is-info}\n\n', string)
    pattern = re.compile(r'(\w+)::\s*\n(.*?)\s*\n\s*\n', re.DOTALL)
    string = re.sub(pattern, r'\1:\n```python\n\2\n```\n\n', string)
    pattern = re.compile(r'^(.*?)\s*={3,}\s*$', re.MULTILINE)
    string = re.sub(pattern, r'## \1\n', string)
    string = transform_parm_and_return(string)
    return string

def process_one_file(file):
    path = file.relative_to("..")
    if file.stem == "__init__":
        path = file.parent.relative_to("..")
    name = str(path.with_suffix("")).replace("\\", ".").replace("/", ".")
    m = Module(name, path)
    current_working = m
    current_status = None
    current_parent = m
    current_parent_space = 0
    with open(file, encoding="UTF-8") as f:
        for line in f:
            nspace = len(line) - len(line.lstrip())
            line = line.strip()
            if not current_status and current_working and not line.startswith("'''") and not line.startswith('"""'):
                current_working = None
                current_status = None
            elif not current_status and current_working and line.startswith("'''") and line.endswith("'''") and line != "'''":
                current_working.docstring = line.split("'''")[1].strip()
                current_working = None
                current_status = None
            elif not current_status and current_working and line.startswith('"""') and line.endswith('"""') and line != '"""':
                current_working.docstring = line.split('"""')[1].strip()
                current_working = None
                current_status = None
            elif not current_status and current_working and line.startswith("'''"):
                current_working.docstring = line.strip().split('"""')[1].strip()
                current_status = "'''"
            elif not current_status and current_working and line.startswith('"""'):
                current_working.docstring = line.strip().split('"""')[1].strip()
                current_status = '"""'
            elif current_status and current_working:
                if line.endswith(current_status):
                    current_working.docstring += line.split(current_status)[0]
                    current_working = None
                    current_status = None
                else:
                    current_working.docstring += line + "\n"
            elif not current_working and line.startswith("class"):
                if "(" in line:
                    name = line.split("(")[0].split("class")[1].strip()
                else:
                    name = line.split(":")[0].split("class")[1].strip()
                if name.startswith("_"):
                    continue
                current_working = Class(name)
                m.classes.append(current_working)
                current_parent = current_working
                current_parent_space = nspace
            elif not current_working and line.startswith("def"):
                name = line.split("(")[0].split("def")[1].strip()
                if name.startswith("_"):
                    continue
                current_working = Function(name)
                if nspace == 0 and current_parent == m:
                    current_parent.functions.append(current_working)
                elif nspace > current_parent_space:
                    current_parent.functions.append(current_working)
                else:
                    current_parent = m
                    current_parent_space = 0
    return m

xponge = Path("../Xponge")
paths = {}

modules = []

for root, dirs, files in os.walk(xponge):
    root = Path(root)
    if root.stem == '__pycache__':
        continue
    modules.append(process_one_file(root / "__init__.py"))
    for file in files:
        file = root / file
        if file.stem.startswith("__") or not file.suffix.endswith(".py"):
            continue
        modules.append(process_one_file(file))
        

for root, dirs, files in os.walk(xponge):
    root = Path(root)
    if root.stem == '__pycache__':
        continue
    parent = root.parent
    parent = parent.relative_to("..")
    try:
        parent = parent.with_suffix("")    
        that = str(parent).replace("\\", ".").replace("/", ".")
        this = str(root.relative_to("..").with_suffix("").as_posix())
        if that in Module.modules:
            Module.modules[that].subpackages.append(this)
    except ValueError:
        pass

    that = str(root.relative_to("..").as_posix()).replace("\\", ".").replace("/", ".")
    for file in files:
        file = root / file
        if file.stem.startswith("__") or not file.suffix.endswith(".py"):
            continue
        this = str(file.relative_to("..").with_suffix("").as_posix())
        if that in Module.modules:
            Module.modules[that].submodules.append(this)

for module in modules:
    module_folder = "Xponge文档" / module.path.parent
    module_folder.mkdir(parents=True, exist_ok=True)
    with open("Xponge文档" / module.path.with_suffix(".md"), "w", encoding="utf-8") as f:
        print(f"""# {module.name}

{process_docstring(module.docstring)}""", file=f)
        print("\n## subpackages\n", file=f)
        for subpackage in module.subpackages:
            print(f"\n### [{subpackage.replace('/', '.')}](/文档/Xponge文档/{subpackage})\n", file=f)
        print("\n## submodules\n", file=f)
        for submodule in module.submodules:
            print(f"\n### [{submodule.replace('/', '.')}](/文档/Xponge文档/{submodule})\n", file=f)
        print("\n## functions\n", file=f)
        for function in module.functions:
            if not function.docstring:
                continue
            print(f"""### {function.name}

{process_docstring(function.docstring)}
""", file=f)
        if module.classes:
            print("## classes\n", file=f)
        for class_ in module.classes:
            print(f"""### {class_.name}

{process_docstring(class_.docstring)}
""", file=f)
            for function in class_.functions:
                if not function.docstring:
                    continue
                print(f"""#### {function.name}

{process_docstring(function.docstring)}
""", file=f)
