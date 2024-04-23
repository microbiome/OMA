# Use a specific version of nixpkgs from an bleeding-edge fork of github.com/NixOS/nixpkgs
let
  pkgs = import (fetchTarball "https://github.com/rstats-on-nix/nixpkgs/archive/bc30fa6c937a4598564bd0fa7114cbad2a3554aa.tar.gz") {};
  # Add generic R packages required for this build
  rpkgs = with pkgs.rPackages; [
    BiocManager
    BiocBook
  ];

  # Build mia package
  mia = [(pkgs.rPackages.buildRPackage {
     name = "mia";
     src = pkgs.fetchgit {
      url = "https://github.com/microbiome/mia";
      branchName = "devel";
      rev = "58cc49bd118ff2b75d12b2eca006cfd48c307933";
      sha256 = "sha256-BsQ7/nJBLadF2NQgGpHhAC10AFQaiK+HbMCIjfOH6d4=";
  };
    # mia dependencies (see DESCRIPTION)
    propagatedBuildInputs = builtins.attrValues {
    inherit (pkgs.rPackages)
      SummarizedExperiment
      SingleCellExperiment
      TreeSummarizedExperiment
      MultiAssayExperiment
      MASS
      ape
      decontam
      vegan
      BiocGenerics
      S4Vectors
      IRanges
      Biostrings
      DECIPHER
      BiocParallel
      DelayedArray
      DelayedMatrixStats
      scuttle
      scater
      DirichletMultinomial
      rlang
      dplyr
      tibble
      tidyr
      bluster
      rgl
      MatrixGenerics;
    };
  }) ];

  # Build miaTime package
  miatime = [(pkgs.rPackages.buildRPackage {
      name = "miaTime";
      src = pkgs.fetchgit {
        url = "https://github.com/microbiome/miaTime";
        branchName = "master";
        rev = "9fe9771f7329fc991796eb79cc1e17ee06e1bc24";
        sha256 = "sha256-IL9CbL0HWKlpmMHq1Rxen7+utzpw4qfb4NuVvM0N0oA=";
    };
    # miaTime dependencies (see DESCRIPTION)
    propagatedBuildInputs = builtins.attrValues {
      inherit (pkgs.rPackages)
      dplyr
      S4Vectors
      SummarizedExperiment
      SingleCellExperiment
      vegan;
    } ++ [ mia ];
  }) ];

  # Build SpiecEasi package
  spieceasi = [(pkgs.rPackages.buildRPackage {
     name = "SpiecEasi";
     src = pkgs.fetchgit {
      url = "https://github.com/zdk123/SpiecEasi";
      branchName = "master";
      rev = "5f396da85baa114b31c13d9744c05387a1b04c23";
      sha256 = "sha256-Z3x7hK2ieLxjQVn94DCPJCDP86TK+k5no4/e5jb8ihg=";
  };
    # SpiecEasi dependencies (see DESCRIPTION)
    propagatedBuildInputs = builtins.attrValues {
    inherit (pkgs.rPackages)
      huge
      pulsar
      MASS
      VGAM
      Matrix
      glmnet
      RcppArmadillo;
    };
  }) ];

  # Build SPRING package
  spring = [(pkgs.rPackages.buildRPackage {
     name = "SPRING";
     src = pkgs.fetchgit {
      url = "https://github.com/GraceYoon/SPRING";
      branchName = "master";
      rev = "3d641a4b939b1b3cc042c064a05000aa48266af0";
      sha256 = "sha256-H1kEy5dPPjiUPFiQLFzbdsO5t204NSCPnQfqPQitMTs=";
  };
    # SPRING dependencies (see DESCRIPTION)
    propagatedBuildInputs = builtins.attrValues {
    inherit (pkgs.rPackages)
      mixedCCA
      huge
      pulsar
      rootSolve
      mvtnorm;
    } ++ [spieceasi];
  }) ];

  # Build NetCoMi package
  netcomi = [(pkgs.rPackages.buildRPackage {
     name = "NetCoMi";
     src = pkgs.fetchgit {
      url = "https://github.com/stefpeschel/NetCoMi";
      branchName = "main";
      rev = "0809c7a5e0f1e74cb9023fbf1186d477739cc6f7";
      sha256 = "sha256-X+isckPsojo2rfaICXXmydN1AcT1IzOJaMqGEe9CIxE=";
  };
    # NetCoMi dependencies (see DESCRIPTION)
    propagatedBuildInputs = builtins.attrValues {
    inherit (pkgs.rPackages)
      Biobase
      corrplot
      doSNOW
      fdrtool
      filematrix
      foreach
      gtools
      huge
      igraph
      MASS
      Matrix
      mixedCCA
      orca
      phyloseq
      pulsar
      qgraph
      RColorBrewer
      Rdpack
      rlang
      vegan
      WGCNA;
    } ++ [spring];
  }) ];

  # Build OMA book/package
  oma = [(pkgs.rPackages.buildRPackage {
     name = "oma";
     src = pkgs.fetchgit {
       url = "https://github.com/microbiome/OMA";
       branchName = "devel";
       rev = "67dd77ef36f7b90c416ef98ce5c8f2086f64fbdb";
       sha256 = "sha256-U4fHmVKvoCq6YNeSgoT1eY1ZD6AiX+xiAXtq9pQl5Ak=";
     };
     # oma dependencies
     propagatedBuildInputs = builtins.attrValues {
       inherit (pkgs.rPackages)
        rebook
        glue
        sessioninfo
        microbiomeDataSets
        microbiome
        ggsignif
        SummarizedExperiment
        TreeSummarizedExperiment
        kableExtra
        dendextend
        NbClust
        randomcoloR
        cobiclust
        biclust
        tidyverse
        ALDEx2
        ANCOMBC
        Maaslin2
        MicrobiomeStat
        GUniFrac
        devtools
        ComplexHeatmap
        mikropml
        MLeval
        sechm
        ggpubr
        fido
        rgl
        miaViz;
     } ++ [ miatime spieceasi];
   })
];

   # System dependencies
   system_packages = builtins.attrValues {
    inherit (pkgs) R glibcLocales quarto ;
 };

  # R wrapper for nix
  R = pkgs.rWrapper.override {
    packages = rpkgs ++ miatime ++ mia ++ spieceasi ++ spring ++ netcomi ++ oma;
  };

in
# Build R environment
pkgs.mkShell {
  LOCALE_ARCHIVE = if pkgs.system == "x86_64-linux" then  "${pkgs.glibcLocalesUtf8}/lib/locale/locale-archive" else "";
  LANG = "en_US.UTF-8";
  LC_ALL = "en_US.UTF-8";
  LC_TIME = "en_US.UTF-8";
  LC_MONETARY = "en_US.UTF-8";
  LC_PAPER = "en_US.UTF-8";
  LC_MEASUREMENT = "en_US.UTF-8";

  buildInputs = [
    R
    system_packages
  ];
}
