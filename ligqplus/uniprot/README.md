uniprot\_tools
==============

uniprot provides a command-line and python interface to access the
uniprot database

available services: map, retrieve

**map**
    map a list of ids from one format onto another using uniprots mapping api
    
    Args:
        query: id or list of ids to be mapped
        f: from ACC | P_ENTREZGENEID | ...
        t: to ...
        format: tab by default

    Help:
        for a list of all possible mappings visit
        'https://www.uniprot.org/help/api_idmapping'
    

**retrieve**
    request entries by uniprot acc using batch retrieval

    Args:
        query: list of ids to retrieve
        format: fasta by default

    Help:
        possible formats:
        txt, xml, rdf, gff

Installation
------------

### From pypi (recommended)

    pip install uniprot_tools

### From source (UNIX) as standalone only

Clone the git repository
    
    git clone https://github.com/fedeserral/uniprot.git
 ##Modification of:
https://github.com/jdrudolph/uniprot.git

Use `distutils` to install the package

    cd uniprot
    sudo python setup.py install

Example
-------

### standalone

    uniprot map ACC P_ENTREZGENEID acc_file map_file

This will read UniprotIDs seperated by whitespaces from `acc_file` and
store them to `map_file`.

    uniprot retrieve acc_file entries.txt

Retrieve textual etries for all uniprot ids in `acc_file` and save to
`entries.txt`

Using a pipe:

    echo P31749 | uniprot map ACC P_ENTREZGENEID

will print the result to `stdout` which can be redirected further

    echo P31749 | uniprot retrieve

will print the result to `stdout` which can be redirected further

### inside a python script

    import uniprot as uni
    print uni.map('P31749', f='ACC', t='P_ENTREZGENEID') # map single id
    print uni.map(['P31749','Q16204'], f='ACC', t='P_ENTREZGENEID') # map list of ids
    print uni.retrieve('P31749')
    print uni.retrieve(['P31749','Q16204'])
