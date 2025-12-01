from modeller import *
from modeller.automodel import *

# 1. Define Sequences
# --------------------------------------------------
# Lysozyme N-term (Res 1-10 in lys.pdb)
lys_n_seq = "MNIFEMLRID"

# Insert from binder (Res 117-132 in binder_only.pdb)
# Identified by matching sequence FQEF... in binder_only.pdb ATOM records [cite: 248-296]
binder_seq = "FQEFSPNLWGLEFQKN"

# Lysozyme C-term (Res 23-164 in lys.pdb)
# Starts with GYYT... which corresponds to Gly 23 in lys.pdb 
lys_c_seq = "".join([
    "GYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCNGVITK",
    "DEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRCALINMVFQMGETGVAGFTNSLRM",
    "LQQKRWDAAAAALAAAAWYNQTPNRAKRVITTFRTGTWDAYKNL",
])

# Full Target Sequence
target_seq = lys_n_seq + binder_seq + lys_c_seq

# 2. Create Alignment File (.ali)
# --------------------------------------------------
ali_file = open('chimera.ali', 'w')

# --- Entry 1: Chimera Target ---
ali_file.write(">P1;chimera\n")
ali_file.write("sequence:chimera:.:.:.:.::::\n")
ali_file.write("%s*\n" % target_seq)

# --- Entry 2: Template 1 (Lysozyme N-term) ---
# Source: lys.pdb | Chain: A | Range: 1-10
ali_file.write(">P1;template_lys_n\n")
ali_file.write("structureX:lys:1:A:10:A::::\n")
ali_file.write("%s%s%s*\n" % (lys_n_seq, '-' * len(binder_seq), '-' * len(lys_c_seq)))

# --- Entry 3: Template 2 (Binder Loop) ---
# Source: binder_only.pdb | Chain: B | Range: 117-132
ali_file.write(">P1;template_binder\n")
ali_file.write("structureX:binder_only:117:B:132:B::::\n")
ali_file.write("%s%s%s*\n" % ('-' * len(lys_n_seq), binder_seq, '-' * len(lys_c_seq)))

# --- Entry 4: Template 3 (Lysozyme C-term) ---
# Source: lys.pdb | Chain: A | Range: 23-164
ali_file.write(">P1;template_lys_c\n")
ali_file.write("structureX:lys:23:A:164:A::::\n")
ali_file.write("%s%s%s*\n" % ('-' * len(lys_n_seq), '-' * len(binder_seq), lys_c_seq))

ali_file.close()

# 3. Run Automodel
# --------------------------------------------------
log.verbose()
env = environ()
env.io.atom_files_directory = ['.']

a = automodel(env,
              alnfile  = 'chimera.ali',
              knowns   = ('template_lys_n', 'template_binder', 'template_lys_c'),
              sequence = 'chimera')

a.starting_model = 1
a.ending_model   = 50
a.make()