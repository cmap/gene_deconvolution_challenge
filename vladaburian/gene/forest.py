
from gene import gio, conf


experiment = conf.EXPERIMENTS[0]
ds = gio.load_prepared_dataset(conf.RES, experiment)

