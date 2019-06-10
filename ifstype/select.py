def xval(x):
    "Choose the interval containing x"
    def select(iv_net):
        for iv in iv_net.net:
            if x in iv:
                return iv
    return select

def random(iv_net):
    "Choose a random interval"
    return iv_net.random_net_iv()[0]

def min_delta(iv_net):
    "Choose the interval of the smallest size"
    return iv_net.argmin_delta()
