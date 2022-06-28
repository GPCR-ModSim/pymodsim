class Printing(object):
    def __init__(self):
        """
        This is a proxy to put a message directly to the stdout through
        *print* command
        """
        pass

    def dispatch(self, msg):
        """
        Simply print the msg passed
        """
        print(msg)