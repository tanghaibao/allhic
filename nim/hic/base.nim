import future
import logging
import strutils
import "matrix"


var logger: Logger
let loggers = getHandlers()

##
## Make sure there is only one instance of logger
##
if loggers.len >= 1:
    logger = loggers[0]
else:
    logger = newConsoleLogger()
    logger.addHandler()