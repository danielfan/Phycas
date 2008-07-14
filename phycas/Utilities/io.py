import os
_test_data_dir = None
def getPhycasTestDataDir():
    global _test_data_dir
    if _test_data_dir is None:
        d = os.path.dirname(os.path.dirname(__file__))
        _test_data_dir = os.path.join(d, "Tests", "Data")
    return _test_data_dir

def getPhycasTestData(filen):
    return os.path.join(getPhycasTestDataDir(), filen)

