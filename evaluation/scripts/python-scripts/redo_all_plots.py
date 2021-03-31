import os
import sys
import re

def exec():
    plots_path = sys.argv[0][:-len(sys.argv[0].split('/')[-1])] + "../../plots/"
    filepaths = os.listdir(plots_path)
    tests = {}
    for path in filepaths:
        if os.path.isdir(plots_path + path): tests[path] = {}
    for test in tests:
        test_scenes = []
        for scene in os.listdir(plots_path + test):
            if not os.path.isdir(plots_path + test + '/' + scene): continue
            if not re.match(r"((crown|measure-one|villa|hair|landscape|ecosys):)*(crown|measure-one|villa|hair|landscape|ecosys)", scene):
                continue
            test_scenes.append(scene.replace(':', ','))
        tests[test] = test_scenes
    for test in tests.copy():
        if not tests[test]: del tests[test]
    rappy_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(rappy_path)
    with open("../redo_all_plots.sh", 'w') as f:
        f.write("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" &> /dev/null && pwd )\"\n")
        f.write("cd $DIR/../results\n")
        for test in tests:
            f.write("cd " + test + '\n')
            for scene in tests[test]:
                command = "./run.sh " + scene + " --skip-render"
                f.write(command + '\n')
            f.write("cd ..\n")
    os.system("./../redo_all_plots.sh")

exec()