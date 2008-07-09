import sys, os
import textwrap
import phycas.ReadNexus as ReadNexus
from phycas import getDefaultOutFilter, OutFilter

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


class PhycasHelp(object):
    _phycas_cmd_classes = set()
    def __str__(self):
        PhycasTablePrinter._reset_term_width()
        n = [i.__name__.lower() for i in PhycasHelp._phycas_cmd_classes]
        n.sort()
        done = False
        i = iter(n)
        o = []
        while not done:
            a, b, c = ("", "", "")
            try:
                a = i.next()
                b = i.next()
                c = i.next()
            except StopIteration:
                done = True
            r = PhycasTablePrinter.format_help(a, b, c)
            o.append(r)
        t = "\n".join(o)
        return """Phycas Help
    
For Python Help use "python_help()"

  Commands objects are called by following the name by () and then
   hitting the RETURN key. Thus, to invoke the sumt command use:

sumt()

  Commands (and almost everything else in python) are case-sensitive -- so
    "Sumt" is _not_ the same thing as "sumt"

  Lower case versions are the instances of command objects that most users 
    will be interested in using.

  The currently implemented Phycas commands are:
  
%s

  Use <command_name>.help to see the detailed help for each command. So,

sumt.help

  will display the help information for the sumt command object.
""" % t

    def __call__(self):
        """Returns the result of __str__, so that users can omit the () on the
        and simply type "cmd.help" to see the help message."""
        print(str(self))

    def __repr__(self):
        return str(self)
phycas_help = PhycasHelp()

class FileFormats:
    NEXUS, FASTA, PHYLIP, PDF, RAW_TEXT = range(5)
    _names = ["NEXUS", "FASTA", "PHYLIP", "PDF", "RAW_TEXT"]
    def to_str(v):
        if v < 0 or v > len(FileFormats._names):
            raise VauleError("Invalid file format code (%s) specified" % str(v))
        return FileFormats._names[v]
    to_str = staticmethod(to_str)

class ExistingFileBehavior(object):
    pass

class ReplaceExistingFileBehavior(ExistingFileBehavior):
    def __str__(self):
        return "REPLACE"

class AppendExistingFileBehavior(ExistingFileBehavior):
    def __str__(self):
        return "APPEND"

class AddNumberExistingFileBehavior(ExistingFileBehavior):
    def __str__(self):
        return "ADD_NUMBER"

REPLACE = ReplaceExistingFileBehavior()
APPEND = AppendExistingFileBehavior()
ADD_NUMBER = AddNumberExistingFileBehavior()

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

_opt_name_help_len = 24
_opt_val_help_len = 19

class PhycasTablePrinter:
    _help_str_wrapper = textwrap.TextWrapper()
    _text_wrapper = textwrap.TextWrapper()
    _help_fmt_str = "%%-%ds %%-%ds %%s" % (_opt_name_help_len, _opt_val_help_len)

    def _reset_term_width():
        tty_sz = ttysize()
        if tty_sz:
            PhycasTablePrinter._set_terminal_width(tty_sz[1])
    _reset_term_width = staticmethod(_reset_term_width)

    def _set_terminal_width(w):
        global _opt_val_help_len, _opt_name_help_len
        x = _opt_val_help_len + _opt_name_help_len + 2
        new_width = max(x + _opt_name_help_len + 1, w)
        PhycasTablePrinter._help_str_wrapper.width = new_width
        PhycasTablePrinter._help_str_wrapper.subsequent_indent = " "*x
        PhycasTablePrinter._help_str_wrapper.initial_indent = " "*x
        PhycasTablePrinter._text_wrapper.width = new_width
        PhycasTablePrinter._text_wrapper.subsequent_indent = "    "
        PhycasTablePrinter._text_wrapper.initial_indent = ""
    _set_terminal_width = staticmethod(_set_terminal_width)

    def format_help(name, value, helpStr):
        global _opt_val_help_len, _opt_name_help_len
        hs = "\n".join(PhycasTablePrinter._help_str_wrapper.wrap(helpStr))
        w = _opt_val_help_len + _opt_name_help_len + 2
        return PhycasTablePrinter._help_fmt_str % (name, value, hs[w:])
    format_help = staticmethod(format_help)

    def get_help_divider():
        global _opt_val_help_len, _opt_name_help_len
        l = PhycasTablePrinter._help_str_wrapper.width - (_opt_val_help_len + _opt_name_help_len + 2)
        return " ".join(["="*_opt_name_help_len, "="*_opt_val_help_len, "="*l])
    get_help_divider = staticmethod(get_help_divider)
PhycasTablePrinter._set_terminal_width(80)

class PhycasOutput(object):
    def __init__(self):
        self._is_active = True
    def _silence(self):
        self._is_active = False
    def _activate(self):
        self._is_active = True

class FileOutputSpec(PhycasOutput):
    def __init__(self, prefix="", help_str="", filename=None):
        PhycasOutput.__init__(self)
        self.mode = ADD_NUMBER
        self.filename = filename
        self.prefix = prefix
        if filename and prefix:
            self.prefix = None
        self._help_str = help_str
        self._valid_formats = []
        self._options = []
        self._opened_filename = None

    def _help_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called 
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = ""
        v = None
        fn = self._getFilename()
        if fn:
            if self._is_active:
                v = fn
            else:
                v = "<silenced>"
        opts_help = [PhycasTablePrinter.format_help("%s" %pref, v, self._help_str)]
        if fn:
            if self.prefix:
                opts_help.append(PhycasTablePrinter.format_help("%sprefix" % dpref, self.prefix, "file prefix (appropriate suffix will be added)"))
            if self.filename:
                opts_help.append(PhycasTablePrinter.format_help("%sfilename" % dpref, self.filename, "The full file name."))
            opts_help.append(PhycasTablePrinter.format_help("%smode" % dpref, str(self.mode), 'Controls the behavior when the file is present. Valid settings are %s. ADD_NUMBER indicates that a number will be added to the end of the file name (or prefix) to make the name unique' % self._getValidModeNames()))
            if len(self._valid_formats) > 1:
                opts_help.append(PhycasTablePrinter.format_help("%s.format" % dpref, FileFormats.to_str(self.format), 'Format of the file valid settings are %s' % self._getValidFormatNames()))
        if self._options:
            opts_help.extend(self._options._help_str_list(pref))
        return opts_help

    def __setattr__(self, name, v):
        if name.startswith("_"):
            self.__dict__[name] = v
        elif name == "filename":
            if v:
                self.__dict__["prefix"] = None
            self.__dict__["filename"] = v
        elif name == "prefix":
            if v:
                self.__dict__["filename"] = None
            self.__dict__["prefix"] = v
        elif name == "mode":
            if not isinstance(v, ExistingFileBehavior):
                raise ValueError('Valid settings of the "mode" attribute are %s'  % self._getValidModeNames())
            self.__dict__["mode"] = v
        elif name == "format":
            if not (v  in self._valid_formats):
                raise ValueError('Valid settings of the "format" attribute are %s' % self._getValidFormatNames())
            self.__dict__["format"] = v
        elif self._options and name in self._options:
            self._options._set_opt(name.lower(), v)
        else:
            raise AttributeError('Output specifier has no "%s" attribute"' % name)

    def _getValidFormatNames(self):
        v = self._valid_formats
        n = [FileFormats.to_str(i) for i in v]
        lv = len(v)
        if lv == 0:
            return None
        if lv == 1:
            return n[0]
        if lv == 2:
            return "%s or %s" % (n[0], n[1])
        return "%s, or %s" % (", ".join(n[0:-1]), n[-1])

    def set(self, v):
        s = str(v)
        if s == v:
            self.filename = s
        else:
            raise ValueError("An output setting can (currently) only be set to a string that is taken to be the filename")
    def _getValidModeNames(self):

        return 'REPLACE, APPEND, or ADD_NUMBER'
    def __str__(self):
       PhycasTablePrinter._reset_term_width()
       return "\n".join(self._help_str_list())

    def _getFilename(self):
        if self._opened_filename:
            return self._opened_filename
        pref, suffix = None, None
        if self.filename:
            fn = self.filename
            pref, suffix = os.path.splitext(fn)
        elif self.prefix:
            pref, suffix = self.prefix, self._getSuffix()
        if pref:
            fn = "%s%s" % (pref, suffix)
            if isinstance(self.mode, AddNumberExistingFileBehavior):
                curr_index = 1
                while os.path.exists(fn):
                    fn = "%s%d%s" % (pref, curr_index ,suffix)
                    curr_index += 1
            return fn
        return None

class TextOutputSpec(FileOutputSpec):
    def __init__(self, prefix="", help_str="", filename=None):
        FileOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.RAW_TEXT]

class BinaryOutputSpec(FileOutputSpec):
    def __init__(self, prefix="", help_str="", filename=None):
        FileOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.RAW_TEXT]
    def _getValidModeNames(self):
        return 'REPLACE, or ADD_NUMBER'

class TreeOutputSpec(TextOutputSpec):
    def __init__(self, prefix="", help_str="", filename=None):
        FileOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.NEXUS]
        self.format = FileFormats.NEXUS

    def _getSuffix(self):
        return ".tre"

    
class PDFOutputSpec(BinaryOutputSpec):
    def _getSuffix(self):
        return ".pdf"
    def __init__(self, prefix="", help_str="", filename=None):
        BinaryOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.PDF]
        self.format = FileFormats.PDF

class PhycasCommandOutputOptions(object):
    """The PhycasCommandOutputOptions is a very simple class written to
    provide a standardized place for command output settings to reside so 
    that users can clearly see what the output options are (without cluttering
    command settings with lots of settings such as tree_out_replace, 
    tree_out_append...
    """
    def __init__(self, verbosity_level=None):
        self.__dict__["level"] = 0
        self.__dict__["_cached_level"] = 0
        if verbosity_level is None:
            self.level = getDefaultOutFilter()
        else:
            self.level = verbosity_level
        self.__dict__["_help_order"] = []
    def __setattr__(self, name, value):
        o = self.__dict__.get(name)
        if o is not None:
            if name == "level":
                if value < OutFilter.DEBUGGING or value > OutFilter.SILENT:
                    raise ValueError("The verbosity level must be set to one of the following:\n    OutFilter.%s" % "\n    OutFilter.".join(["DEBUGGING", "VERBOSE", "NORMAL", "WARNINGS", "ERRORS", "SILENT"]))
                self.__dict__["level"] = value
                self.__dict__["_cached_level"] = value
            else:
                isintarg = isinstance(value, int) or isinstance(value, long)
                turning_off = ((value is None) or (value is False) or (isintarg and value == 0))
                turning_on = ((value is True) or (isintarg and value != 0))
                if turning_off and (outp is not None):
                    o._silence()
                elif turning_on and (outp is not None):
                    o._activate()
                else:
                    o.set(value)
        else:
            if name in OutFilter._names:
                i = OutFilter._names.index(name)
                if value:
                    if self.level > value:
                        self.__dict__["level"] = value
                        self.__dict__["_cached_level"] = value
                elif self.level < value + 1:
                    self.__dict__["level"] = value
                    self.__dict__["_cached_level"] = value
            else:
                raise AttributeError("%s has no attribute %s" % (self.__class__.__name__, name))

    def _silence(self):
        """Silence the output from the associated function.
        
        Somewhat counterintuitively, this function moves the level down
        to OutFilter.ERRORS rather than OutFilter.SILENT.
        The OutFilter.SILENT level is intended only for developers 
        It can be set by (for an instance called `out`):
            out._silence()
            out.level = OutFilter.SILENT
        """
        self.__dict__["level"] = OutFilter.ERRORS
        for n in self.__dict__["_help_order"]:
            self.__dict__[n]._silence()
       
    def _activate(self):
        self.__dict__["level"] = self._cached_level
        for n in self.__dict__["_help_order"]:
            self.__dict__[n]._silence()

    def _help_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called 
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = pref
        s = "OutFilter." + OutFilter.to_str(self.level)
        opts_help = [PhycasTablePrinter.format_help("%slevel" %dpref, s, "Controls the amount of output (verbosity) of the command")]
        for n in self.__dict__["_help_order"]:
            a = self.__dict__[n]
            opts_help.append("")
            opts_help.extend(a._help_str_list("%s%s" % (dpref, n)))
        return opts_help

    def __str__(self):
       PhycasTablePrinter._reset_term_width()
       return "\n".join(self._help_str_list())



class PhycasCmdOpts(object):
    
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
    def _help_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called 
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = ""
        opts_help = []
        for i in self._optionsInOrder:
            n = pref and "%s%s" % (dpref, i[0]) or i[0]
            s = PhycasTablePrinter.format_help(n, self._current[i[0]], i[2])
            opts_help.append(s)
        return opts_help

    def __str__(self):
        PhycasTablePrinter._reset_term_width()
        return "\n".join(self._help_str_list())

    def _set_opt(self, name, value):
        if name in self._current:
            transf = self._transformer.get(name)
            if transf is None:
                self._current[name] = value
            else:
                try:
                    self._current[name] = transf(self, value)
                except:
                    try:
                        c = ""
                        print transf
                        c = transf.get_description()
                        print c
                    except:
                        raise
                    raise ValueError("%s is not a valid value for %s.  %s" % (value, name, c))
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
    def get_description(self):
        if self.min is None:
            if self.max is None:
                return "Must be an integer."
            else:
                return "Must be an integer no greater than %d" % self.max
        else:
            if self.max is None:
                return "Must be an integer no less than %d" % self.min
            else:
                return "Must be an integer in the range from %d to %d " % (self.min, self.max)

class FloatArgValidate(object):
    def __init__(self, min=None, max=None):
        self.min = min
        self.max = max
    def __call__(self, opts, v):
        iv = float(v)
        if (self.min is not None) and (self.min > iv):
            raise ValueError("Value must be >= %s" % str(self.min))
        if (self.max is not None) and (self.max < iv):
            raise ValueError("Value must be <= %s" % str(self.min))
        return iv
    def get_description(self):
        if self.min is None:
            if self.max is None:
                return "Must be a real number."
            else:
                return "Must be a real number no greater than %f" % self.max
        else:
            if self.max is None:
                return "Must be a real number no less than %f" % self.min
            else:
                return "Must be a real number in the range from %f to %f " % (self.min, self.max)

class ProbArgValidate(FloatArgValidate):
    def __init__(self):
        FloatArgValidate.__init__(self, min=0.0, max=1.0)

class PhycasCommandHelp(object):
    def __init__(self, command, cmd_name, cmd_descrip):
        self.command = command
        self.cmd_name = cmd_name
        self.cmd_descrip = cmd_descrip

    def __str__(self):
        PhycasTablePrinter._reset_term_width()
        command = self.command
        opts = command.__dict__["_options"]
        out = command.__dict__["out"]
        d = PhycasTablePrinter.get_help_divider()
        a = [self.cmd_name]
        a.extend(PhycasTablePrinter._text_wrapper.wrap(self.cmd_descrip))
        o = opts._help_str_list()
        h = PhycasTablePrinter.format_help("Attribute", "Current Value", "Explanation")
        if len(o) > 0:
            a.extend(["\nAvailable options:", h, d])
            a.extend(o)
            a.extend([d, ""])
        if out is not None:
            a.extend(["Current output settings:", h, d])
            a.extend(out._help_str_list("out"))
        return "\n".join(a)

    def __call__(self):
        """Returns the result of __str__, so that users can omit the () on the
        and simply type "cmd.help" to see the help message."""
        print(str(self))

    def __repr__(self):
        return str(self)

class PhycasCommand(object):

    def __init__(self, phycas, option_defs, output_options=None):
        # the roundabout way of initializing PhycasCmdOpts is needed because
        #   _options must be in the __dict__ before the setattr is called
        #   on a PhycasCommand
        PhycasHelp._phycas_cmd_classes.add(self.__class__)
        o = PhycasCmdOpts()
        self.__dict__["_options"] = o
        o._initialize(self, option_defs)
        self.__dict__["phycas"] = phycas
        self.__dict__["out"] = output_options
        self.__dict__["help"] = PhycasCommandHelp(self, "Sumt", "The Sumt command is used to summarize a collection of trees (usually trees that have been produced by an MCMC simulation).")

    def __setattr__(self, name, value):
        o = self.__dict__["_options"]
        if name in o:
            o._set_opt(name.lower(), value)
        elif name in self.__dict__:
            outp = self.__dict__["out"]
            isintarg = isinstance(value, int) or isinstance(value, long)
            turning_off = ((value is None) or (value is False) or (isintarg and value == 0))
            turning_on = ((value is True) or (isintarg and value != 0))
            if name == "out":
                if turning_off and (outp is not None):
                    outp._silence()
                elif turning_on and (outp is not None):
                    outp._activate()
                else:
                    self.__dict__["out"] = value
            else:
                self.__dict__[name] = value
        else:
            raise AttributeError("%s has no attribute %s" % (self.__class__.__name__, name))

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
