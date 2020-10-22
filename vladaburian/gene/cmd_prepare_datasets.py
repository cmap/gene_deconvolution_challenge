from gene import gio, conf


for experiment in conf.EXPERIMENTS:
    gio.cmd_prepare_dataset(conf.RES, experiment)

print('DONE')
