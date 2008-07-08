import sys, os
import textwrap
import phycas.ReadNexus as ReadNexus

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

class FileFormats:
    NEXUS, FASTA, PHYLIP = range(3)

class TreeSource(object):
    def __init__(self, arg, **kwargs):
        """`arg` can be a string or an iterable containing trees.
        If `arg` is a string, then the `format` keyword argument can be used to
        specify the file format (the default is NEXUS).
        If the trees are passed in as a list then a `taxon_labels` keyword
        argument is expected.
        """
        if isinstance(arg, str):
            self.filename = arg
            self.format = kwargs.get("format", FileFormats.NEXUS)
            self.trees = None
            self.taxon_labels = None
        else:
            self.trees = list(arg)
            self.taxon_labels = kwargs.get("taxon_labels")
            self.filename = None
            self.format = None
            
    def description(arg):
        if isinstance(arg, TreeSource) and arg.filename:
            return 'trees from the file "%s"' % arg.filename
        return "collection of trees"
    description = staticmethod(description)

    def __iter__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns an iterable over the trees (this will trigger the reading of 
        the input file if the object was initialized from a string.
        Stores tree descriptions in self.stored_tree_defs and taxon labels in
        self.taxon_labels.
        """
        if self.trees is None:
            if not self.filename:
                raise ValueError("No file name specified for the source of trees")
            if not os.path.exists(self.filename):
                raise ValueError('The file "%s" does not exist')
            reader = ReadNexus.NexusReader()
            reader.readFile(self.filename)
            self.taxon_labels = reader.getTaxLabels()
            self.trees = reader.getTrees()
        return iter(self.trees)


class PhycasCmdOpts(object):
    _help_str_wrapper = textwrap.TextWrapper()
    _help_fmt_str = "%-19s %-19s %s"
    
    def __init__(self, command=None, args=None):
        self.__dict__["_current"] = {}
        self.__dict__["_transformer"] = {}
        self.__dict__["_help_info"] = {}
        self.__dict__["_optionsInOrder"] = None
        self.__dict__["_command"] = None
        self.__dict__["_unchecked"] = set()
        if command and args:
            self._initialize(command, args)
    def _initialize(self, command=None, args=None):
        self._current = {}
        self._transformer = {}
        self._help_info = {}
        self._optionsInOrder = args
        self._command = command
        for opt in args:
            if len(opt) == 3:
                name, default, help_str = opt
                transf = None
            else:
                name, default, help_str, transf = opt
            self._help_info[name] =  [default, help_str]
            self._current[name] = default
            if transf is not None:
                self._transformer[name] = transf
            self._command.__dict__[name] = self._current[name]
    def __str__(self):
        PhycasCmdOpts._reset_term_width()
        opts_help = []
        for i in self._optionsInOrder:
            hs = "\n".join(PhycasCmdOpts._help_str_wrapper.wrap(i[2]))
            s = PhycasCmdOpts._help_fmt_str % (i[0], self._current[i[0]], hs[40:])
            opts_help.append(s)
        return "\n".join(opts_help)
    def _set_opt(self, name, value):
        if name in self._current:
            transf = self._transformer.get(name)
            if transf is None:
                self._current[name] = value
            else:
                try:
                    self._current[name] = transf(self, value)
                except:
                    raise ValueError("%s is not a valid value for %s" % (value, name))
            self._command.__dict__[name] = self._current[name]
        else:
            raise AttributeError("%s does not contain an attribute %s" % (self._command.__class__.__name__, name))
    def __contains__(self, name):
        return name in self._current
    def __getattr__(self, name):
        if name in self._current:
            return self._current[name]
        raise AttributeError("%s does not contain an attribute %s" % (self.__class__.__name__, name))
    def __setattr__(self, name, value):
        if name in self.__dict__:
            self.__dict__[name] = value
        else:
            self._set_opt(name, value)
    def check_unchecked(self):
        for name in self._unchecked:
            value = self._current[name]
            try:
                transf = self._transformer.get(name)
                if not transf is None:
                    self._current[name] = transf(self, value)
            except:
                pass # it is concievable that some errors here are not errors once all of the other opts have been transformed, so we set options twice
        for name in self._unchecked:
            value = self._current[name]
            self._set_opt(name, value)
        self._unchecked.clear()
    def set_unchecked(self, name, value):
        "Sets the attribute `key` to the `value`"
        if name in self._current:
            self._current[name] = value
            self._unchecked.add(name)
        else:
            raise AttributeError("%s does not contain an attribute %s" % (self._command.__class__.__name__, name))
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


def TreeSourceValidate(opts, v):
    return TreeSource(v)

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

class PhycasCommandHelp(object):
    def __init__(self, cmd_name, cmd_descrip, opts):
        self.opts = opts
        self.cmd_name = cmd_name
        self.cmd_descrip = cmd_descrip
    def __str__(self):
        return "\n".join([self.cmd_name, self.cmd_descrip, "Current settings:", str(self.opts)])
    def __call__(self):
        print(str(self))
    def __repr__(self):
        return str(self)

class PhycasCommand(object):
    def __init__(self, phycas, option_defs):
        # the roundabout way of initializing PhycasCmdOpts is needed because
        #   _options must be in the __dict__ before the setattr is called
        #   on a PhycasCommand
        o = PhycasCmdOpts()
        self.__dict__["_options"] = o
        o._initialize(self, option_defs)
        self.__dict__["phycas"] = phycas
        self.__dict__["help"] = PhycasCommandHelp("Sumt", "The Sumt command is used to summarize a collection of trees (usually trees that have been produced by an MCMC simulation).", o)
    def __getattr__(self, name):
        raise AttributeError("%s does not contain an attribute %s" % (self.__class__.__name__, name))
    def __setattr__(self, name, value):
        o = self.__dict__["_options"]
        if name in o:
            o._set_opt(name.lower(), value)
        elif name in self.__dict__:
            self.__dict__[name] = value
        else:
            raise AttributeError("%s does not contain an attribute %s" % (self.__class__.__name__, name))

    def set(self, **kwargs):
        "Sets the attributes of a command using keyword arguments"
        old = {}
        o = self._options
        try:
            for key, value in kwargs.iteritems():
                old[key] = getattr(self, key)
                o.set_unchecked(key, value)
            o.check_unchecked()
        except:
            for key, value in old.iteritems():
                o.set_unchecked(key, value)
            raise
