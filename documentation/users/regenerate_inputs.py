import re

categories = {
    'mcmc': ('random_seed','ncycles','sample_every','report_every','outfile_prefix','data_source','data_file_name','starting_tree_source','tree_topology','nchains'),
    #'sumt': ('sumt_outgroup_taxon','trees','sumt_trees_prefix','sumt_splits_prefix','sumt_output_replace','sumt_burnin','sumt_equal_brlens'),
    'edges': ('internal_edgelen_dist','external_edgelen_dist','edgelen_dist','using_hyperprior','edgelen_hyperprior')
    #'xxxx': ('','','','','','','','','','','','',''),
    }

re_init = re.compile(r'def __init__\(self\):(.*)def source_line\(\):', re.MULTILINE | re.DOTALL)
re_entry = re.compile(r'self\.(\S+)\s+=\s+(\S+)\s+# (.*)')
re_underscore = re.compile('_')
re_bracketed = re.compile(r"<(\S+)>")
re_quoted = re.compile(r"'(\S+)'")
re_dblquoted = re.compile(r'"(\S+)"')
re_true = re.compile('True')
re_false = re.compile('False')
re_none = re.compile('None')

def italtt(m):
    return "{\\tt \\em %s}" % m.group(1)

def justdblquoted(m):
    return "``%s''" % m.group(1)

def quotedsmalltt(m):
    return "{\\tt \\small '%s'}" % m.group(1)

def justsmalltt(m):
    return "{\\tt \\small %s}" % m.group(0)

if __name__ == '__main__':
    contents = open('../../phycas/Phycas/Phycas.py','r').read()
    terms = {}
    match = re_init.search(contents)
    if match:
        lines = match.group(1).split('\n')
        for line in lines:
            m = re_entry.match(line.strip())
            if m:
                varname = m.group(1)
                vardefinition = m.group(3)
                vardefault = m.group(2)
                terms[varname] = (vardefault,vardefinition)
    else:
        print 'sorry, no match'

    term_names = terms.keys()
    term_names.sort()

    # create all.tex containing an alphabetical list of all settings
    #outf = open('all.tex','w')
    #outf.write('\\begin{description}\n')
    #for name in term_names:
    #    default, definition = terms[name]
    #    varname = re_underscore.sub(r'\_', name)
    #    vardefinition = re_underscore.sub(r'\_', definition)
    #    vardefault = re_underscore.sub(r'\_', default)
    #    outf.write('\\item[%s] %s (default: {\\tt %s})\n' % (varname, vardefinition, vardefault))
    #outf.write('\\end{description}\n')

    for cat in categories.keys():
        print 'category',cat,':',categories[cat]
        # create a category-specific alphabetical list of settings
        outf = open('%s.tex' % cat,'w')
        outf.write('\\begin{description}\n')
        for name in term_names:
            if name in categories[cat]:
                default, definition = terms[name]
                varname = re_underscore.sub(r'\_', name)
                vardefinition = re_underscore.sub(r'\_', definition)
                vardefinition = re_quoted.sub(quotedsmalltt, vardefinition)
                vardefinition = re_dblquoted.sub(justdblquoted, vardefinition)
                vardefinition = re_true.sub(justsmalltt, vardefinition)
                vardefinition = re_false.sub(justsmalltt, vardefinition)
                vardefinition = re_none.sub(justsmalltt, vardefinition)
                vardefinition = re_bracketed.sub(italtt, vardefinition)
                vardefault = re_underscore.sub(r'\_', default)
                outf.write('\\item[%s] \\index{%s} %s (default: {\\tt \\small %s})\n' % (varname, varname, vardefinition, vardefault))
        outf.write('\\end{description}\n')
