# Localisation du SDK Gmsh ; pour l'évaluation il sera bien dans le dossier parent
GMSH_DIR := ../../../../../gmsh-sdk-linux

# Choix du compilateur
CC := gcc

# Flags de compilation ; vous pouvez ajouter d'autres flags si vous le souhaitez
CFLAGS := -Wall

# Chemins vers les dossiers `include` où se trouvent les headers externes
# on donne le chemin vers le dossier `include` du SDK
INC_DIR := -I $(GMSH_DIR)/include -I header # specify include directories

# Chemins vers les dossiers `lib` où se trouvent les librairies externes
# on donne le chemin vers le dossier `lib` du SDK
LIB_FILE := $(GMSH_DIR)/lib/libgmsh.so.4.11.1 # specify library file for compilation

# Spécification du runtime path
LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib -lm

# On indique au compilateur de linker la librairie Gmsh
LDLIBS := -lgmsh # -llapack

# Nom du programme à compiler
PROG := project

# Liste des objets nécessaires pour compiler le programme
# à modifier si vous ajoutez d'autres modules !
OBJS := elasticity.o lu.o matrix.o design.o project.o eigen.o

# Règle de compilation
all: $(PROG)

# Règle de compilation : on produit un .o à partir d'un .c
# on fournit les flags de compilation CFLAGS et les chemins include INC_DIR
# `$<` signifie: "la première dépendance" (le .c)
%.o: %.c
	$(CC) -g -c $(CFLAGS) $(INC_DIR) $<

# Règle de link : on produit le programme PROG à partir des fichiers objets OBJS
# on fournit la librairie à linker et le runtime path
# `$^` signifie : "toutes les dépendances" (les OBJS)
$(PROG): $(OBJS)
	$(CC) -g -o $@ $^ $(LIB_FILE) $(LDFLAGS)

# Règle de nettoyage : supprime PROG et OBJS
clean:
	rm $(PROG) $(OBJS)

savior:
	gcc -o project elasticity.o lu.o matrix.o design.o project.o eigen.o ../../../gmsh-sdk/lib/gmsh-4.11.dll -W -lm