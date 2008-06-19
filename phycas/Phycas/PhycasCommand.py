import sys
import textwrap


###############################################################################
#   Taken from
#   http://mail.python.org/pipermail/python-list/2006-February/365594.html

def ttysize():
    try:
        import fcntl
        import struct
        import termios
        buf = 'abcdefgh'
        buf = fcntl.ioctl(0, termios.TIOCGWINSZ, buf)
        row, col, rpx, cpx = struct.unpack('hhhh', buf)
        return row, col
    except:
        return None
###############################################################################

class PhycasCmdOpts(object):
    _help_str_wrapper = textwrap.TextWrapper()
    _help_fmt_str = "%-19s %-19s %s"
    
    def __init__(self, command, args):
        self.current = {}
        self.transformer = {}
        self.help_info = {}
        self.optionsInOrder = args
        self.command = command
        for opt in args:
            if len(opt) == 3:
                name, default, help_str = opt
                transf = None
            else:
                name, default, help_str, transf = opt
            self.help_info[name] =  [default, help_str]
            self.current[name] = default
            if transf is not None:
                self.transformer[name] = transf

    def __str__(self):
        PhycasCmdOpts._reset_term_width()
        opts_help = []
        for i in self.optionsInOrder:
            hs = "\n".join(PhycasCmdOpts._help_str_wrapper.wrap(i[2]))
            s = PhycasCmdOpts._help_fmt_str % (i[0], self.current[i[0]], hs[40:])
            opts_help.append(s)
        return "\n".join(opts_help)
    def set_opt(self, name, value):
        if name in self.current:
            transf = self.transformer.get(name)
            if transf is None:
                self.current[name] = value
            else:
                try:
                    self.current[name] = transf(self, value)
                except:
                    raise #ValueError("%s is not a valid value for %s" % (value, name))
        else:
            raise AttributeError("%s does not contain an attribute %s" % (self.command.__class__.__name__, name))
    def _reset_term_width():
        tty_sz = ttysize()
        if tty_sz:
            PhycasCmdOpts._set_terminal_width(tty_sz[1])
    _reset_term_width = staticmethod(_reset_term_width)

    def _set_terminal_width(w):
        new_width = max(60, w)
        PhycasCmdOpts._help_str_wrapper.width = new_width
        PhycasCmdOpts._help_str_wrapper.subsequent_indent = " "*40
        PhycasCmdOpts._help_str_wrapper.initial_indent = " "*40
    _set_terminal_width = staticmethod(_set_terminal_width)

PhycasCmdOpts._set_terminal_width(80)

def BoolArgValidate(opts, v):
    return bool(v)
class IntArgValidate(object):
    def __init__(self, min=None, max=None):
        self.min = min
        self.max = max
    def __call__(self, opts, v):
        iv = int(v)
        if (self.min is not None) and (self.min > iv):
            raise ValueError("Value must be >= %s" % str(self.min))
        if (self.max is not None) and (self.max < iv):
            raise ValueError("Value must be <= %s" % str(self.min))
        return iv

class PhycasCommand():
    def __init__(self, phycas, option_defs):
        self.__dict__["options"] = PhycasCmdOpts(self, option_defs)
        self.__dict__["phycas"] = phycas
    def __getattr__(self, name):
        try:
            return self.options.current[name.lower()]
        except:
            raise AttributeError("%s does not contain an attribute %s" % (self.__class__.__name__, name))
    def __setattr__(self, name, value):
        if name in self.__dict__:
            self.__dict__[name] = value
        else:
            self.__dict__["options"].set_opt(name.lower(), value)
    def help(self, opt=None):
       print self.options
        
