import timeit

startTime = None
prevTime = None
prevStr = None

def log(str, printStartTime=False):
    global startTime, prevTime, prevStr
    if startTime is None:
        startTime = timeit.default_timer()
        prevTime= timeit.default_timer()

    t = timeit.default_timer()

    fromStart = t - startTime
    fromPrev = t - prevTime

    if prevStr is not None:
        print(f'{prevStr} finished - {fromPrev: 1.3f} s')

    if printStartTime:
        print(f'{str} started at {fromStart: 1.3f} s')
    else:
        print(f'{str} started')

    prevTime = t
    prevStr = str
