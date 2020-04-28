#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os.path

class Params:
    def __init__(self, params=None, scriptpath=sys.argv[0], withargument=None, required=None, allowed=None, helpoption=None, usage=None, help=None):
        if None == params:
          params = sys.argv[1:]
        if None == scriptpath:
          scriptpath = sys.argv[0]
        
        self.scriptpath = scriptpath
        self.scriptname = os.path.basename(scriptpath)
        self.scriptdir  = os.path.dirname(scriptpath)
        
        self.usage = usage
        self.help  = help
        self.helpoption = helpoption
        self.parse(params, withargument, required, allowed)
    
    # Make object seem like a list
    #
    def __iter__(self):
        return self.data.__iter__()
    def __getitem__(self, i):
        return self.data.__getitem__(i)
    def __getslice__(self, i, j):
        return self.data.__getslice__(i,j)
    def __len__(self):
        return self.data.__len__()
    def __str__(self):
        return self.data.__str__()
    
    
    def parse(self, data, withargument, required, allowed):
        self.data = data
        self.args = list(data)
        self.opts = {}
        
        if withargument or required:
          if not allowed:
            allowed = set([])
          else:
            allowed = set(allowed)
          
          if withargument:
            allowed.update(withargument)
          if required:
            allowed.update(required)
        self.withargument = withargument
        self.allowed = allowed
        self.required = required
        
        
        while self.args and self.args[0][0]=='-':
            if self.args[0] == '-':
                break
            
            arg = self.args.pop(0)
            if arg[1]=='-':
                # Force end of option parsing?
                #
                if (arg == "--"):
                    break
                
                # Long options
                #
                arg = arg[2:].split('=', 1)
                opt = arg[0]
                optarg = None
                if len(arg)>1:
                  optarg = arg[1]
                elif None != withargument and opt in withargument:
                  try:
                    optarg = self.args.pop(0)
                  except IndexError:
                    sys.stderr.write("Argument error: option '%s' requires an argument.\n" % (opt))
                    self.writeUsage()
                    sys.exit(1)
                self.opts[opt] = optarg
            else:
                # Short options
                #
                prev = None
                for i in range(1,len(arg)):
                    opt = arg[i]
                    
                    if opt == '=':
                      if not prev:
                        sys.stderr.write("Argument error: invalid option: '='.\n")
                        self.writeUsage()
                        sys.exit(1)
                      optarg = arg[i+1:]
                      self.opts[prev] = optarg
                      break
                      
                    elif None != withargument and opt in withargument:
                      if i < len(arg)-1:
                        optarg = arg[i+1:]
                      else:
                        try:
                          optarg = self.args.pop(0)
                        except IndexError:
                          sys.stderr.write("Argument error: option '%s' requires an argument.\n" % (opt))
                          self.writeUsage()
                          sys.exit(1)
                      self.opts[opt] = optarg
                      break
                    
                    self.opts[opt] = None
                    prev = opt
        
        if self.helpoption in self.opts:
          self.writeHelp()
          sys.exit(0)
        
        self.checkAllowedOptions()
        self.checkOptionArguments()
        self.checkRequiredOptions()
        
    
    def checkAllowedOptions(self):
        if None != self.allowed:
          for opt in self.opts:
            if opt not in self.allowed:
              sys.stderr.write("Argument error: unknown option '%s'.\n" % (opt))
              self.writeUsage()
              sys.exit(1)
    
    
    def checkRequiredOptions(self):
        if None != self.required:
          for opt in self.required:
            if opt not in self.opts:
              sys.stderr.write("Argument error: missing option '%s'.\n" % (opt))
              self.writeUsage()
              sys.exit(1)
    
    
    def checkOptionArguments(self):
        if None != self.withargument:
          for opt in self.opts:
            arg = self.opts[opt]
            if opt in self.withargument:
              if None == arg:
                sys.stderr.write("Argument error: option '%s' requires an argument.\n" % (opt))
                self.writeUsage()
                sys.exit(1)
            else:
              if None != arg:
                sys.stderr.write("Argument error: option '%s' does not take an argument.\n" % (opt))
                self.writeUsage()
                sys.exit(1)
    
    
    def _format_option_usage(self, o):
      if len(o) > 1:
        dashes = "--"
      else:
        dashes = "-"
      
      if self.withargument and o in self.withargument:
        arg = " <%s_value>"%o
      else:
        arg = ""
      
      return dashes + o + arg
    
    
    def writeUsage(self, outstream = sys.stderr, minimal=False):
        if None == self.usage:
          req = ""
          if self.required:
            for o in sorted(self.required):
              req += " " + self._format_option_usage(o)
          
          options = set([])
          if self.allowed:
            options.update(self.allowed)
          if self.required:
            options.difference_update(self.required)
          
          outstream.write("USAGE:\n\t%s [OPTIONS]%s <ARGUMENTS>\n" % (self.scriptname, req))
          outstream.write("\nOPTIONS:\n")
          for o in sorted(options):
            outstream.write("\t%s\n"%(self._format_option_usage(o)))
          
        else:
          outstream.write(self.usage.replace("$0", self.scriptname))
          if self.usage[-1:] != "\n":
            outstream.write("\n")
        
        if not minimal and self.helpoption:
          if len(self.helpoption) > 1:
            dashes = "--"
          elif len(self.helpoption) == 1:
            dashes = "-"
          else:
            dashes = ""
          outstream.write("\nFOR MORE HELP, TYPE:\n\t%s %s%s\n"%(self.scriptname, dashes, self.helpoption))
    
    
    def writeHelp(self, outstream = sys.stdout):
        self.writeUsage(outstream, minimal=True)
        if self.help:
          outstream.write(self.help.replace("$0", self.scriptname))
          if self.help[-1:] != "\n":
            outstream.write("\n")
        
    
    def doUpperCaseDisablesOption(self):
        options = list(self.opts.keys())
        for opt in options:
            if not opt.isupper():
                continue
            
            del self.opts[opt]
            
            opt_lower = c.lower()
            if opt_lower in self.opts:
                del self.opts[opt_lower]
    
    
    def setDefaultOptions(self, **defaults):
        for k in defaults:
          if k not in self.opts:
            self.opts[k] = defaults[k]
    
    
    def getOpt(self, name, default=None):
        if name in self.opts:
          return self.opts[name]
        return default
    
    
    def isOpt(self, name):
        return name in self.opts
    
    
    def getArgs(self, first, num, defaults=[], formats=[]):
        offset = num - len(defaults)
        assert offset >= 0, "Cannot provide more default values than requested arguments!"
        
        output = []
        val = None
        for i in range(num):
          
          if i+first < len(self.args):
            # Get command line argument (string)
            val = self.args[i+first]
          else:
            # Get default value
            assert i >= offset, "No argument at position %d and no default value provided." % (i)
            val = defaults[i-offset]
          
          # Format the value (e.g. type conversion)
          if formats and len(formats) > i and formats[i]:
            f = formats[i]
            call = f.__call__ # checks if f is a callable object, e.g. a function
            val = f(val)
          
          #Output the value
          output.append(val)
        
        return output
    
    
