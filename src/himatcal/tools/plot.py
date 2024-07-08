
def import_style():
    import matplotlib.pyplot as plt
    import pkg_resources

    try:
        mplstyle_file = pkg_resources.resource_filename("himatcal", "tools/science-1.mplstyle")
        plt.style.use(mplstyle_file)
        plt.rcParams['font.family'] = 'Calibri, Microsoft YaHei'
        print("Matplotlib style loaded!")
    except:
        print("font not found")
        pass
