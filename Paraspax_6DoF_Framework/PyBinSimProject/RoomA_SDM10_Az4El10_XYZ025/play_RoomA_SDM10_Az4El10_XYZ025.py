import pybinsim
import logging

pybinsim.logger.setLevel(logging.WARNING)    # defaults to INFO
#Use logging.WARNING for printing warnings only

with pybinsim.BinSim('config/RoomA_SDM10.cfg') as binsim:
    binsim.stream_start()
