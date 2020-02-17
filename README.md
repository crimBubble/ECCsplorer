# ECCsplorer

## Trouble shooting

### segemehl v0.3.4

htslib/sam.h not found:
install htslib from https://github.com/samtools/htslib/releases/tag/1.10.2

in segemehl-0.3.4/include/segemehl.h
change line 54 to: #include "/path_to/htslib/include/htslib/sam.h"

htslib.pc not found:
copy htslib.pc from /path_to/htslib/lib/pkgconfig to /usr/lib/pkgconfig
