#
# Copyright John Reid 2007
#

"""
Code to help with logging and process priority
"""

_logging_setup = False
def setup_process(root_dir, make_nice=True):
    import logging, os.path
    global _logging_setup
    logger = logging.getLogger()
    if not _logging_setup:
        if not os.access(root_dir, os.R_OK):
            os.makedirs(root_dir)
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.StreamHandler())
        logger.addHandler(logging.FileHandler(os.path.join(root_dir, 'log.txt'), 'w'))

        if make_nice:
            try:
                import cookbook
                cookbook.make_current_process_nice()
            except:
                logger.warn('Could not set process priority')

        _logging_setup = True

    return logger
