      MODULE CSMVersion
        TYPE VersionType
          INTEGER :: Major = 4
          INTEGER :: Minor = 8
          INTEGER :: Model = 5
          INTEGER :: Build = 0
        END TYPE VersionType
        TYPE (VersionType) Version
        CHARACTER(len=*), PARAMETER :: VBranch = '-release'
      END MODULE CSMVersion
