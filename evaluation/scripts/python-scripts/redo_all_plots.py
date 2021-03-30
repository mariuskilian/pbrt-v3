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
    cwd = os.path.dirname(os.path.realpath(__file__))
    os.system("cd " + cwd)
    os.system("cd ../../results")
    for test in tests:
        os.system("cd " + test)
        for scene in tests[test]:
            os.system("./run.sh " + scene + " --skip-render")
        os.system("cd ..")

exec()