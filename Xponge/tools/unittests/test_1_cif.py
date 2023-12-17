"""
    This **module** includes unittests of the cif files
"""
from io import StringIO

__all__ = ["test_cof",
           "test_som"]

def test_cof():
    """
        Test loading a cif file for Covalent Organic Framework
    """
    import Xponge
    from Xponge.forcefield.amber import gaff
    from Xponge.mdrun import run

    cif = StringIO("""data_358
_audit_creation_date              2023-11-14
_audit_creation_method            'Materials Studio'
_symmetry_space_group_name_H-M    'P1'
_symmetry_Int_Tables_number       1
_symmetry_cell_setting            triclinic
loop_
_symmetry_equiv_pos_as_xyz
  x,y,z
_cell_length_a                    32.9466
_cell_length_b                    57.0651
_cell_length_c                    6.2388
_cell_angle_alpha                 90.0000
_cell_angle_beta                  90.0000
_cell_angle_gamma                 90.0000
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy
C1     C     0.45858   0.76381   0.75000   0.00000  Uiso   1.00
C2     C     0.37882   0.79040   0.75000   0.00000  Uiso   1.00
C3     C     0.12106   0.87632   0.75000   0.00000  Uiso   1.00
C4     C     0.04139   0.90288   0.75000   0.00000  Uiso   1.00
C5     C     0.45304   0.78894   0.75000   0.00000  Uiso   1.00
C6     C     0.41619   0.80150   0.75000   0.00000  Uiso   1.00
C7     C     0.30122   0.79415   0.75000   0.00000  Uiso   1.00
C8     C     0.27068   0.81253   0.75000   0.00000  Uiso   1.00
C9     C     0.27170   0.85515   0.75000   0.00000  Uiso   1.00
C10    C     0.22915   0.85422   0.75000   0.00000  Uiso   1.00
C11    C     0.19862   0.87259   0.75000   0.00000  Uiso   1.00
C12    C     0.11906   0.90054   0.75000   0.00000  Uiso   1.00
C13    C     0.08181   0.91266   0.75000   0.00000  Uiso   1.00
C14    C     0.46509   0.70706   0.75000   0.00000  Uiso   1.00
C15    C     0.46467   0.68236   0.75000   0.00000  Uiso   1.00
C16    C     0.53317   0.62855   0.75000   0.00000  Uiso   1.00
C17    C     0.52088   0.60409   0.75000   0.00000  Uiso   1.00
C18    C     0.45644   0.58328   0.75000   0.00000  Uiso   1.00
C19    C     0.47911   0.56248   0.75000   0.00000  Uiso   1.00
C20    C     0.46682   0.53803   0.75000   0.00000  Uiso   1.00
C21    C     0.46467   0.48427   0.75000   0.00000  Uiso   1.00
C22    C     0.46512   0.45959   0.75000   0.00000  Uiso   1.00
C23    C     0.58188   0.75403   0.75000   0.00000  Uiso   1.00
C24    C     0.61914   0.76617   0.75000   0.00000  Uiso   1.00
C25    C     0.66560   0.82733   0.75000   0.00000  Uiso   1.00
C26    C     0.70844   0.83341   0.75000   0.00000  Uiso   1.00
C27    C     0.77187   0.81159   0.75000   0.00000  Uiso   1.00
C28    C     0.79174   0.83333   0.75000   0.00000  Uiso   1.00
C29    C     0.83456   0.83941   0.75000   0.00000  Uiso   1.00
C30    C     0.91627   0.86521   0.75000   0.00000  Uiso   1.00
C31    C     0.95307   0.87778   0.75000   0.00000  Uiso   1.00
O32    O     0.29284   0.77339   0.75000   0.00000  Uiso   1.00
O33    O     0.20700   0.89335   0.75000   0.00000  Uiso   1.00
O34    O     0.56850   0.63474   0.75000   0.00000  Uiso   1.00
O35    O     0.43149   0.53183   0.75000   0.00000  Uiso   1.00
O36    O     0.63866   0.84190   0.75000   0.00000  Uiso   1.00
O37    O     0.86152   0.82484   0.75000   0.00000  Uiso   1.00
N38    N     0.34016   0.80329   0.75000   0.00000  Uiso   1.00
N39    N     0.15970   0.86344   0.75000   0.00000  Uiso   1.00
H40    H     0.47274   0.80244   0.75000   0.00000  Uiso   1.00
H41    H     0.41810   0.82035   0.75000   0.00000  Uiso   1.00
H42    H     0.28811   0.87157   0.75000   0.00000  Uiso   1.00
H43    H     0.14636   0.91094   0.75000   0.00000  Uiso   1.00
H44    H     0.09223   0.92925   0.75000   0.00000  Uiso   1.00
H45    H     0.43498   0.71016   0.75000   0.00000  Uiso   1.00
H46    H     0.43544   0.67388   0.75000   0.00000  Uiso   1.00
H47    H     0.42359   0.58328   0.75000   0.00000  Uiso   1.00
H48    H     0.43542   0.49272   0.75000   0.00000  Uiso   1.00
H49    H     0.93502   0.95650   0.75000   0.00000  Uiso   1.00
H50    H     0.59228   0.73742   0.75000   0.00000  Uiso   1.00
H51    H     0.64647   0.75579   0.75000   0.00000  Uiso   1.00
H52    H     0.78829   0.79517   0.75000   0.00000  Uiso   1.00
H53    H     0.91822   0.84636   0.75000   0.00000  Uiso   1.00
H54    H     0.47276   0.36427   0.75000   0.00000  Uiso   1.00
C55    C     0.95858   0.26381   0.75000   0.00000  Uiso   1.00
C56    C     0.87882   0.29040   0.75000   0.00000  Uiso   1.00
C57    C     0.62106   0.37632   0.75000   0.00000  Uiso   1.00
C58    C     0.54139   0.40288   0.75000   0.00000  Uiso   1.00
C59    C     0.95304   0.28894   0.75000   0.00000  Uiso   1.00
C60    C     0.91619   0.30150   0.75000   0.00000  Uiso   1.00
C61    C     0.80122   0.29415   0.75000   0.00000  Uiso   1.00
C62    C     0.77068   0.31253   0.75000   0.00000  Uiso   1.00
C63    C     0.77170   0.35515   0.75000   0.00000  Uiso   1.00
C64    C     0.72915   0.35422   0.75000   0.00000  Uiso   1.00
C65    C     0.69862   0.37259   0.75000   0.00000  Uiso   1.00
C66    C     0.61906   0.40054   0.75000   0.00000  Uiso   1.00
C67    C     0.58181   0.41266   0.75000   0.00000  Uiso   1.00
C68    C     0.96509   0.20706   0.75000   0.00000  Uiso   1.00
C69    C     0.96467   0.18236   0.75000   0.00000  Uiso   1.00
C70    C     0.03317   0.12855   0.75000   0.00000  Uiso   1.00
C71    C     0.02088   0.10409   0.75000   0.00000  Uiso   1.00
C72    C     0.95644   0.08328   0.75000   0.00000  Uiso   1.00
C73    C     0.97911   0.06248   0.75000   0.00000  Uiso   1.00
C74    C     0.96682   0.03803   0.75000   0.00000  Uiso   1.00
C75    C     0.96467   0.98427   0.75000   0.00000  Uiso   1.00
C76    C     0.96512   0.95959   0.75000   0.00000  Uiso   1.00
C77    C     0.08188   0.25403   0.75000   0.00000  Uiso   1.00
C78    C     0.11914   0.26617   0.75000   0.00000  Uiso   1.00
C79    C     0.16560   0.32733   0.75000   0.00000  Uiso   1.00
C80    C     0.20844   0.33341   0.75000   0.00000  Uiso   1.00
C81    C     0.27187   0.31159   0.75000   0.00000  Uiso   1.00
C82    C     0.29174   0.33333   0.75000   0.00000  Uiso   1.00
C83    C     0.33456   0.33941   0.75000   0.00000  Uiso   1.00
C84    C     0.41627   0.36521   0.75000   0.00000  Uiso   1.00
C85    C     0.45307   0.37778   0.75000   0.00000  Uiso   1.00
O86    O     0.79284   0.27339   0.75000   0.00000  Uiso   1.00
O87    O     0.70700   0.39335   0.75000   0.00000  Uiso   1.00
O88    O     0.06850   0.13474   0.75000   0.00000  Uiso   1.00
O89    O     0.93149   0.03183   0.75000   0.00000  Uiso   1.00
O90    O     0.13866   0.34190   0.75000   0.00000  Uiso   1.00
O91    O     0.36152   0.32484   0.75000   0.00000  Uiso   1.00
N92    N     0.84016   0.30329   0.75000   0.00000  Uiso   1.00
N93    N     0.65970   0.36344   0.75000   0.00000  Uiso   1.00
H94    H     0.97274   0.30244   0.75000   0.00000  Uiso   1.00
H95    H     0.91810   0.32035   0.75000   0.00000  Uiso   1.00
H96    H     0.78811   0.37157   0.75000   0.00000  Uiso   1.00
H97    H     0.64636   0.41094   0.75000   0.00000  Uiso   1.00
H98    H     0.59223   0.42925   0.75000   0.00000  Uiso   1.00
H99    H     0.93498   0.21016   0.75000   0.00000  Uiso   1.00
H100   H     0.93544   0.17388   0.75000   0.00000  Uiso   1.00
H101   H     0.92359   0.08328   0.75000   0.00000  Uiso   1.00
H102   H     0.93542   0.99272   0.75000   0.00000  Uiso   1.00
H103   H     0.43502   0.45650   0.75000   0.00000  Uiso   1.00
H104   H     0.09228   0.23742   0.75000   0.00000  Uiso   1.00
H105   H     0.14647   0.25579   0.75000   0.00000  Uiso   1.00
H106   H     0.28829   0.29517   0.75000   0.00000  Uiso   1.00
H107   H     0.41822   0.34636   0.75000   0.00000  Uiso   1.00
H108   H     0.97276   0.86427   0.75000   0.00000  Uiso   1.00
C217   C     0.54142   0.76381   0.75000   0.00000  Uiso   1.00
C218   C     0.62118   0.79040   0.75000   0.00000  Uiso   1.00
C219   C     0.87894   0.87632   0.75000   0.00000  Uiso   1.00
C220   C     0.95861   0.90288   0.75000   0.00000  Uiso   1.00
C221   C     0.54696   0.78894   0.75000   0.00000  Uiso   1.00
C222   C     0.58381   0.80150   0.75000   0.00000  Uiso   1.00
C223   C     0.69878   0.79415   0.75000   0.00000  Uiso   1.00
C224   C     0.72932   0.81253   0.75000   0.00000  Uiso   1.00
C225   C     0.72830   0.85515   0.75000   0.00000  Uiso   1.00
C226   C     0.77085   0.85422   0.75000   0.00000  Uiso   1.00
C227   C     0.80138   0.87259   0.75000   0.00000  Uiso   1.00
C228   C     0.88094   0.90054   0.75000   0.00000  Uiso   1.00
C229   C     0.91819   0.91266   0.75000   0.00000  Uiso   1.00
C230   C     0.53491   0.70706   0.75000   0.00000  Uiso   1.00
C231   C     0.53533   0.68236   0.75000   0.00000  Uiso   1.00
C232   C     0.46683   0.62855   0.75000   0.00000  Uiso   1.00
C233   C     0.47912   0.60409   0.75000   0.00000  Uiso   1.00
C234   C     0.54356   0.58328   0.75000   0.00000  Uiso   1.00
C235   C     0.52089   0.56248   0.75000   0.00000  Uiso   1.00
C236   C     0.53318   0.53803   0.75000   0.00000  Uiso   1.00
C237   C     0.53533   0.48427   0.75000   0.00000  Uiso   1.00
C238   C     0.53488   0.45959   0.75000   0.00000  Uiso   1.00
C239   C     0.41812   0.75403   0.75000   0.00000  Uiso   1.00
C240   C     0.38086   0.76617   0.75000   0.00000  Uiso   1.00
C241   C     0.33440   0.82733   0.75000   0.00000  Uiso   1.00
C242   C     0.29156   0.83341   0.75000   0.00000  Uiso   1.00
C243   C     0.22813   0.81159   0.75000   0.00000  Uiso   1.00
C244   C     0.20826   0.83333   0.75000   0.00000  Uiso   1.00
C245   C     0.16544   0.83941   0.75000   0.00000  Uiso   1.00
C246   C     0.08373   0.86521   0.75000   0.00000  Uiso   1.00
C247   C     0.04693   0.87778   0.75000   0.00000  Uiso   1.00
O248   O     0.70716   0.77339   0.75000   0.00000  Uiso   1.00
O249   O     0.79300   0.89335   0.75000   0.00000  Uiso   1.00
O250   O     0.43150   0.63474   0.75000   0.00000  Uiso   1.00
O251   O     0.56851   0.53183   0.75000   0.00000  Uiso   1.00
O252   O     0.36134   0.84190   0.75000   0.00000  Uiso   1.00
O253   O     0.13848   0.82484   0.75000   0.00000  Uiso   1.00
N254   N     0.65984   0.80329   0.75000   0.00000  Uiso   1.00
N255   N     0.84030   0.86344   0.75000   0.00000  Uiso   1.00
H256   H     0.52726   0.80244   0.75000   0.00000  Uiso   1.00
H257   H     0.58190   0.82035   0.75000   0.00000  Uiso   1.00
H258   H     0.71189   0.87157   0.75000   0.00000  Uiso   1.00
H259   H     0.85364   0.91094   0.75000   0.00000  Uiso   1.00
H260   H     0.90777   0.92925   0.75000   0.00000  Uiso   1.00
H261   H     0.56502   0.71016   0.75000   0.00000  Uiso   1.00
H262   H     0.56456   0.67388   0.75000   0.00000  Uiso   1.00
H263   H     0.57641   0.58328   0.75000   0.00000  Uiso   1.00
H264   H     0.56458   0.49272   0.75000   0.00000  Uiso   1.00
H265   H     0.06498   0.95650   0.75000   0.00000  Uiso   1.00
H266   H     0.40772   0.73742   0.75000   0.00000  Uiso   1.00
H267   H     0.35353   0.75579   0.75000   0.00000  Uiso   1.00
H268   H     0.21171   0.79517   0.75000   0.00000  Uiso   1.00
H269   H     0.08178   0.84636   0.75000   0.00000  Uiso   1.00
H270   H     0.52724   0.36427   0.75000   0.00000  Uiso   1.00
C271   C     0.04142   0.26381   0.75000   0.00000  Uiso   1.00
C272   C     0.12118   0.29040   0.75000   0.00000  Uiso   1.00
C273   C     0.37894   0.37632   0.75000   0.00000  Uiso   1.00
C274   C     0.45861   0.40288   0.75000   0.00000  Uiso   1.00
C275   C     0.04696   0.28894   0.75000   0.00000  Uiso   1.00
C276   C     0.08381   0.30150   0.75000   0.00000  Uiso   1.00
C277   C     0.19878   0.29415   0.75000   0.00000  Uiso   1.00
C278   C     0.22932   0.31253   0.75000   0.00000  Uiso   1.00
C279   C     0.22830   0.35515   0.75000   0.00000  Uiso   1.00
C280   C     0.27085   0.35422   0.75000   0.00000  Uiso   1.00
C281   C     0.30138   0.37259   0.75000   0.00000  Uiso   1.00
C282   C     0.38094   0.40054   0.75000   0.00000  Uiso   1.00
C283   C     0.41819   0.41266   0.75000   0.00000  Uiso   1.00
C284   C     0.03491   0.20706   0.75000   0.00000  Uiso   1.00
C285   C     0.03533   0.18236   0.75000   0.00000  Uiso   1.00
C286   C     0.96683   0.12855   0.75000   0.00000  Uiso   1.00
C287   C     0.97912   0.10409   0.75000   0.00000  Uiso   1.00
C288   C     0.04356   0.08328   0.75000   0.00000  Uiso   1.00
C289   C     0.02089   0.06248   0.75000   0.00000  Uiso   1.00
C290   C     0.03318   0.03803   0.75000   0.00000  Uiso   1.00
C291   C     0.03533   0.98427   0.75000   0.00000  Uiso   1.00
C292   C     0.03488   0.95959   0.75000   0.00000  Uiso   1.00
C293   C     0.91812   0.25403   0.75000   0.00000  Uiso   1.00
C294   C     0.88086   0.26617   0.75000   0.00000  Uiso   1.00
C295   C     0.83440   0.32733   0.75000   0.00000  Uiso   1.00
C296   C     0.79156   0.33341   0.75000   0.00000  Uiso   1.00
C297   C     0.72813   0.31159   0.75000   0.00000  Uiso   1.00
C298   C     0.70826   0.33333   0.75000   0.00000  Uiso   1.00
C299   C     0.66544   0.33941   0.75000   0.00000  Uiso   1.00
C300   C     0.58373   0.36521   0.75000   0.00000  Uiso   1.00
C301   C     0.54693   0.37778   0.75000   0.00000  Uiso   1.00
O302   O     0.20716   0.27339   0.75000   0.00000  Uiso   1.00
O303   O     0.29300   0.39335   0.75000   0.00000  Uiso   1.00
O304   O     0.93150   0.13474   0.75000   0.00000  Uiso   1.00
O305   O     0.06851   0.03183   0.75000   0.00000  Uiso   1.00
O306   O     0.86134   0.34190   0.75000   0.00000  Uiso   1.00
O307   O     0.63848   0.32484   0.75000   0.00000  Uiso   1.00
N308   N     0.15984   0.30329   0.75000   0.00000  Uiso   1.00
N309   N     0.34030   0.36344   0.75000   0.00000  Uiso   1.00
H310   H     0.02726   0.30244   0.75000   0.00000  Uiso   1.00
H311   H     0.08190   0.32035   0.75000   0.00000  Uiso   1.00
H312   H     0.21189   0.37157   0.75000   0.00000  Uiso   1.00
H313   H     0.35364   0.41094   0.75000   0.00000  Uiso   1.00
H314   H     0.40777   0.42925   0.75000   0.00000  Uiso   1.00
H315   H     0.06502   0.21016   0.75000   0.00000  Uiso   1.00
H316   H     0.06456   0.17388   0.75000   0.00000  Uiso   1.00
H317   H     0.07641   0.08328   0.75000   0.00000  Uiso   1.00
H318   H     0.06458   0.99272   0.75000   0.00000  Uiso   1.00
H319   H     0.56498   0.45650   0.75000   0.00000  Uiso   1.00
H320   H     0.90772   0.23742   0.75000   0.00000  Uiso   1.00
H321   H     0.85353   0.25579   0.75000   0.00000  Uiso   1.00
H322   H     0.71171   0.29517   0.75000   0.00000  Uiso   1.00
H323   H     0.58178   0.34636   0.75000   0.00000  Uiso   1.00
H324   H     0.02724   0.86427   0.75000   0.00000  Uiso   1.00
C433   C     0.50000   0.72240   0.75000   0.00000  Uiso   1.00
C434   C     0.50000   0.66922   0.75000   0.00000  Uiso   1.00
C435   C     0.50000   0.49738   0.75000   0.00000  Uiso   1.00
C436   C     0.50000   0.44427   0.75000   0.00000  Uiso   1.00
N437   N     0.50000   0.64345   0.75000   0.00000  Uiso   1.00
N438   N     0.50000   0.52314   0.75000   0.00000  Uiso   1.00
N439   N     0.50000   0.75001   0.75000   0.00000  Uiso   1.00
C440   C     0.00000   0.22240   0.75000   0.00000  Uiso   1.00
C441   C     0.00000   0.16922   0.75000   0.00000  Uiso   1.00
C442   C     0.00000   0.99738   0.75000   0.00000  Uiso   1.00
C443   C     0.00000   0.94427   0.75000   0.00000  Uiso   1.00
N444   N     0.00000   0.14345   0.75000   0.00000  Uiso   1.00
N445   N     0.00000   0.02314   0.75000   0.00000  Uiso   1.00
N446   N     0.00000   0.25001   0.75000   0.00000  Uiso   1.00
N463   N     0.00000   0.91668   0.75000   0.00000  Uiso   1.00
N464   N     0.50000   0.41668   0.75000   0.00000  Uiso   1.00
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type
C1     C5      1.446   .     D
C1     C239    1.445   .     S
C1     N439    1.576   .     S
C2     C6      1.385   .     D
C2     N38     1.471   .     S
C2     C240    1.384   .     S
C3     C12     1.384   .     D
C3     N39     1.470   .     S
C3     C246    1.384   .     S
C4     C13     1.444   .     D
C4     C247    1.444   .     S
C4     N463    1.575   .     S
C5     C6      1.410   .     S
C5     H40     1.007   .     S
C6     H41     1.078   .     S
C7     C8      1.453   .     S
C7     O32     1.216   .     D
C7     N38     1.385   .     S
C8     C242    1.376   .     S
C8     C243    1.403   .     D
C9     C10     1.403   .     S
C9     H42     1.082   .     S
C9     C242    1.403   .     D
C10    C11     1.453   .     S
C10    C244    1.377   .     D
C11    O33     1.216   .     D
C11    N39     1.385   .     S
C12    C13     1.409   .     S
C12    H43     1.078   .     S
C13    H44     1.007   .     S
C14    C15     1.410   .     D
C14    H45     1.008   .     S
C14    C433    1.445   .     S
C15    H46     1.078   .     S
C15    C434    1.385   .     S
C16    C17     1.453   .     S
C16    O34     1.216   .     D
C16    N437    1.385   .     S
C17    C233    1.376   .     D
C17    C234    1.403   .     S
C18    C19     1.402   .     D
C18    H47     1.082   .     S
C18    C233    1.403   .     S
C19    C20     1.453   .     S
C19    C235    1.377   .     S
C20    O35     1.217   .     D
C20    N438    1.385   .     S
C21    C22     1.408   .     D
C21    H48     1.078   .     S
C21    C435    1.384   .     S
C22    H103    1.007   .     S
C22    C436    1.444   .     S
C23    C24     1.410   .     S
C23    H50     1.008   .     S
C23    C217    1.445   .     D
C24    H51     1.078   .     S
C24    C218    1.384   .     D
C25    C26     1.453   .     S
C25    O36     1.216   .     D
C25    N254    1.385   .     S
C26    C224    1.376   .     S
C26    C225    1.403   .     D
C27    C28     1.403   .     S
C27    H52     1.082   .     S
C27    C224    1.403   .     D
C28    C29     1.453   .     S
C28    C226    1.377   .     D
C29    O37     1.217   .     D
C29    N255    1.384   .     S
C30    C31     1.409   .     S
C30    H53     1.078   .     S
C30    C219    1.384   .     D
C31    H108    1.008   .     S
C31    C220    1.444   .     D
N38    C241    1.385   .     S
N39    C245    1.384   .     S
H49    C76     1.007   .     S
H54    C85     1.008   .     S
C55    C59     1.446   .     D
C55    C293    1.445   .     S
C55    N446    1.576   1_655 D
C56    C60     1.385   .     D
C56    N92     1.471   .     S
C56    C294    1.384   .     S
C57    C66     1.384   .     D
C57    N93     1.470   .     S
C57    C300    1.384   .     S
C58    C67     1.444   .     D
C58    C301    1.444   .     S
C58    N464    1.575   .     S
C59    C60     1.410   .     S
C59    H94     1.007   .     S
C60    H95     1.078   .     S
C61    C62     1.453   .     S
C61    O86     1.216   .     D
C61    N92     1.385   .     S
C62    C296    1.376   .     S
C62    C297    1.403   .     D
C63    C64     1.403   .     S
C63    H96     1.082   .     S
C63    C296    1.403   .     D
C64    C65     1.453   .     S
C64    C298    1.377   .     D
C65    O87     1.216   .     D
C65    N93     1.385   .     S
C66    C67     1.409   .     S
C66    H97     1.078   .     S
C67    H98     1.007   .     S
C68    C69     1.410   .     D
C68    H99     1.008   .     S
C68    C440    1.445   1_655 S
C69    H100    1.078   .     S
C69    C441    1.385   1_655 S
C70    C71     1.453   .     S
C70    O88     1.216   .     D
C70    N444    1.385   .     S
C71    C288    1.403   .     D
C71    C287    1.376   1_455 D
C72    C73     1.402   .     S
C72    H101    1.082   .     S
C72    C287    1.403   .     D
C73    C74     1.453   .     S
C73    C289    1.377   1_655 S
C74    O89     1.217   .     D
C74    N445    1.385   1_655 S
C75    C76     1.408   .     D
C75    H102    1.078   .     S
C75    C442    1.384   1_655 D
C76    C443    1.444   1_655 S
C77    C78     1.410   .     S
C77    H104    1.008   .     S
C77    C271    1.445   .     D
C78    H105    1.078   .     S
C78    C272    1.384   .     D
C79    C80     1.453   .     S
C79    O90     1.216   .     D
C79    N308    1.385   .     S
C80    C278    1.376   .     S
C80    C279    1.403   .     D
C81    C82     1.403   .     S
C81    H106    1.082   .     S
C81    C278    1.403   .     D
C82    C83     1.453   .     S
C82    C280    1.377   .     D
C83    O91     1.217   .     D
C83    N309    1.384   .     S
C84    C85     1.409   .     S
C84    H107    1.078   .     S
C84    C273    1.384   .     D
C85    C274    1.444   .     D
N92    C295    1.385   .     S
N93    C299    1.384   .     S
C217   C221    1.446   .     S
C217   N439    1.576   .     S
C218   C222    1.385   .     S
C218   N254    1.471   .     S
C219   C228    1.384   .     S
C219   N255    1.470   .     S
C220   C229    1.444   .     S
C220   N463    1.575   1_655 S
C221   C222    1.410   .     D
C221   H256    1.007   .     S
C222   H257    1.078   .     S
C223   C224    1.453   .     S
C223   O248    1.216   .     D
C223   N254    1.385   .     S
C225   C226    1.403   .     S
C225   H258    1.082   .     S
C226   C227    1.453   .     S
C227   O249    1.216   .     D
C227   N255    1.385   .     S
C228   C229    1.409   .     D
C228   H259    1.078   .     S
C229   H260    1.007   .     S
C230   C231    1.410   .     S
C230   H261    1.008   .     S
C230   C433    1.445   .     D
C231   H262    1.078   .     S
C231   C434    1.385   .     D
C232   C233    1.453   .     S
C232   O250    1.216   .     D
C232   N437    1.385   .     S
C234   C235    1.402   .     D
C234   H263    1.082   .     S
C235   C236    1.453   .     S
C236   O251    1.217   .     D
C236   N438    1.385   .     S
C237   C238    1.408   .     S
C237   H264    1.078   .     S
C237   C435    1.384   .     D
C238   H319    1.007   .     S
C238   C436    1.444   .     D
C239   C240    1.410   .     D
C239   H266    1.008   .     S
C240   H267    1.078   .     S
C241   C242    1.453   .     S
C241   O252    1.216   .     D
C243   C244    1.403   .     S
C243   H268    1.082   .     S
C244   C245    1.453   .     S
C245   O253    1.217   .     D
C246   C247    1.409   .     D
C246   H269    1.078   .     S
C247   H324    1.008   .     S
H265   C292    1.007   .     S
H270   C301    1.008   .     S
C271   C275    1.446   .     S
C271   N446    1.576   .     S
C272   C276    1.385   .     S
C272   N308    1.471   .     S
C273   C282    1.384   .     S
C273   N309    1.470   .     S
C274   C283    1.444   .     S
C274   N464    1.575   .     S
C275   C276    1.410   .     D
C275   H310    1.007   .     S
C276   H311    1.078   .     S
C277   C278    1.453   .     S
C277   O302    1.216   .     D
C277   N308    1.385   .     S
C279   C280    1.403   .     S
C279   H312    1.082   .     S
C280   C281    1.453   .     S
C281   O303    1.216   .     D
C281   N309    1.385   .     S
C282   C283    1.409   .     D
C282   H313    1.078   .     S
C283   H314    1.007   .     S
C284   C285    1.410   .     S
C284   H315    1.008   .     S
C284   C440    1.445   .     D
C285   H316    1.078   .     S
C285   C441    1.385   .     D
C286   C287    1.453   .     S
C286   O304    1.216   .     D
C286   N444    1.385   1_655 S
C287   C71     1.376   1_655 D
C288   C289    1.402   .     S
C288   H317    1.082   .     S
C289   C290    1.453   .     S
C289   C73     1.377   1_455 S
C290   O305    1.217   .     D
C290   N445    1.385   .     S
C291   C292    1.408   .     S
C291   H318    1.078   .     S
C291   C442    1.384   .     D
C292   C443    1.444   .     D
C293   C294    1.410   .     D
C293   H320    1.008   .     S
C294   H321    1.078   .     S
C295   C296    1.453   .     S
C295   O306    1.216   .     D
C297   C298    1.403   .     S
C297   H322    1.082   .     S
C298   C299    1.453   .     S
C299   O307    1.217   .     D
C300   C301    1.409   .     D
C300   H323    1.078   .     S
C433   N439    1.576   .     S
C434   N437    1.471   .     S
C435   N438    1.470   .     S
C436   N464    1.574   .     S
C440   N446    1.576   .     S
C440   C68     1.445   1_455 S
C441   N444    1.471   .     S
C441   C69     1.385   1_455 S
C442   C75     1.384   1_455 D
C442   N445    1.470   1_565 S
C443   N463    1.574   .     S
C443   C76     1.444   1_455 S
N444   C286    1.385   1_455 S
N445   C74     1.385   1_455 S
N445   C442    1.470   1_545 S
N446   C55     1.576   1_455 D
N463   C220    1.575   1_455 S
""")
    assign, lattice_info = Xponge.get_assignment_from_cif(cif)
    assign.determine_atom_type("gaff")
    assign.calculate_charge("tpacm4")
    restype = assign.to_residuetype("COF")
    gaff.parmchk2_gaff(restype, 'cof.frcmod')
    periodic_bonds = restype.remove_periodic_connectivity()
    lattice = Xponge.Lattice(basis_molecule = restype, periodic_bonds = periodic_bonds, **lattice_info)
    box = Xponge.BlockRegion(0, 0, 0,
        lattice_info["cell_length"][0],
        lattice_info["cell_length"][1],
        100., boundary=True)
    region = Xponge.BlockRegion(0, 0, 0,
        lattice_info["cell_length"][0],
        lattice_info["cell_length"][1],
        lattice_info["cell_length"][2], boundary=True)
    mol = lattice.create(box, region)
    Xponge.save_mol2(mol, "cof.mol2")
    Xponge.save_sponge_input(mol, "cof")
    assert run("SPONGE  -default_in_file_prefix cof -mode NPT -dt 1e-3 \
-thermostat middle_langevin -barostat andersen_barostat -step_limit 5000 > cof1.out") == 0
    box = Xponge.BlockRegion(0, 0, 0,
        2 * lattice_info["cell_length"][0],
        2 * lattice_info["cell_length"][1],
        100., boundary=True)
    region = Xponge.BlockRegion(0, 0, 0,
        2 * lattice_info["cell_length"][0],
        2 * lattice_info["cell_length"][1],
        lattice_info["cell_length"][2], boundary=True)
    mol = lattice.create(box, region)
    Xponge.save_mol2(mol, "cof2.mol2")
    Xponge.save_sponge_input(mol, "cof2")
    assert run("SPONGE  -default_in_file_prefix cof2 -mode NPT -dt 1e-3 \
-thermostat middle_langevin -barostat andersen_barostat -step_limit 5000 > cof2.out") == 0

def test_som():
    """
        Test loading a cif file for Small Organic Molecule
    """
    import Xponge
    from Xponge.forcefield.amber import gaff
    #from Xponge.mdrun import run

    cif = StringIO(r"""data_2100348
loop_
_publ_author_name
'Budzianowski, Armand'
'Katrusiak, Andrzej'
_publ_contact_author_address
;Faculty of Crystal Chemistry,
Adam Mickiewicz University,
Grunwaldzka 6, 60-780 Poznan
Poland
;
_publ_contact_author_email       katran@amu.edu.pl
_publ_contact_author_fax         48618658008
_publ_contact_author_name        'Andrzej Katrusiak'
_publ_contact_author_phone       48618291443
_publ_section_title
;
 Pressure-frozen benzene I revisited
;
_journal_coeditor_code           AV5045
_journal_date_accepted           2005-11-14
_journal_date_recd_electronic    2005-10-03
_journal_issue                   1
_journal_name_full               'Acta Crystallographica Section B'
_journal_page_first              94
_journal_page_last               101
_journal_paper_doi               10.1107/S010876810503747X
_journal_volume                  62
_journal_year                    2006
_chemical_formula_sum            'C6 H6'
_chemical_formula_weight         78.11
_chemical_name_common            benzene
_chemical_name_systematic
;
 benzene
;
_space_group_IT_number           61
_symmetry_cell_setting           orthorhombic
_symmetry_space_group_name_Hall  '-P 2ac 2ab'
_symmetry_space_group_name_H-M   'P b c a'
_atom_sites_solution_hydrogens   geom
_atom_sites_solution_primary     direct
_atom_sites_solution_secondary   difmap
_audit_creation_method           'SHELXL-97 and enCIFer'
_cell_angle_alpha                90.00
_cell_angle_beta                 90.00
_cell_angle_gamma                90.00
_cell_formula_units_Z            4
_cell_length_a                   7.287(6)
_cell_length_b                   9.20(2)
_cell_length_c                   6.688(9)
_cell_measurement_reflns_used    376
_cell_measurement_temperature    296(2)
_cell_measurement_theta_max      29.59
_cell_measurement_theta_min      4.69
_cell_volume                     448.4(12)
_computing_cell_refinement
'CrysAlis RED 1.171.24 beta (Oxford Diffraction Poland, 2004)'
_computing_data_collection
'CrysAlis CCD 1.171.23 beta (Oxford Diffraction Poland, 2004)'
_computing_data_reduction
'CrysAlis RED 1.171.24 beta (Oxford Diffraction Poland, 2004)'
_computing_structure_refinement  'SHELXL-97 (Sheldrick, 1997)'
_computing_structure_solution    'SHELXS-97 (Sheldrick, 1990)'
_diffrn_ambient_temperature      296(2)
_diffrn_detector_area_resol_mean 16.4
_diffrn_measured_fraction_theta_full 0.283
_diffrn_measured_fraction_theta_max 0.283
_diffrn_measurement_device_type  KM4-CCD
_diffrn_measurement_method
;HP omega scans - for more details see:
A. Budzianowski, A. Katrusiak in High-Pressure Crystallography
(Eds.: A. Katrusiak, P. F. McMillan),
Dordrecht: Kluwer Acad. Publ., 2004 pp.157-168
;
_diffrn_radiation_monochromator  graphite
_diffrn_radiation_source         'fine-focus sealed tube'
_diffrn_radiation_type           MoK\a
_diffrn_radiation_wavelength     0.71073
_diffrn_reflns_av_R_equivalents  0.1315
_diffrn_reflns_av_sigmaI/netI    0.1602
_diffrn_reflns_limit_h_max       9
_diffrn_reflns_limit_h_min       -9
_diffrn_reflns_limit_k_max       8
_diffrn_reflns_limit_k_min       -9
_diffrn_reflns_limit_l_max       7
_diffrn_reflns_limit_l_min       -7
_diffrn_reflns_number            1167
_diffrn_reflns_theta_full        29.59
_diffrn_reflns_theta_max         29.59
_diffrn_reflns_theta_min         4.69
_exptl_absorpt_coefficient_mu    0.065
_exptl_absorpt_correction_T_max  0.865
_exptl_absorpt_correction_T_min  0.477
_exptl_absorpt_correction_type   integration
_exptl_absorpt_process_details
;Crystal absorption, DAC absorption and
gasket shadowing absorption has been applied
A. Katrusiak, Z. Kristallogr. 2004, 219, 461-467
;
_exptl_crystal_colour            colourless
_exptl_crystal_density_diffrn    1.157
_exptl_crystal_density_method    'not measured'
_exptl_crystal_F_000             168
_exptl_crystal_size_rad          0.1
_refine_diff_density_max         0.114
_refine_diff_density_min         -0.093
_refine_ls_extinction_method     none
_refine_ls_goodness_of_fit_ref   1.030
_refine_ls_hydrogen_treatment    mixed
_refine_ls_matrix_type           full
_refine_ls_number_parameters     37
_refine_ls_number_reflns         179
_refine_ls_number_restraints     0
_refine_ls_restrained_S_all      1.030
_refine_ls_R_factor_all          0.2023
_refine_ls_R_factor_gt           0.0530
_refine_ls_shift/su_max          0.000
_refine_ls_shift/su_mean         0.000
_refine_ls_structure_factor_coef Fsqd
_refine_ls_weighting_details
'calc w=1/[\s^2^(Fo^2^)+(0.0600P)^2^+0.0300P] where P=(Fo^2^+2Fc^2^)/3'
_refine_ls_weighting_scheme      calc
_refine_ls_wR_factor_gt          0.1097
_refine_ls_wR_factor_ref         0.1355
_reflns_number_gt                94
_reflns_number_total             179
_reflns_threshold_expression     >2sigma(I)
_cod_data_source_file            av5045.cif
_cod_data_source_block           II
_cod_database_code               2100348
loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'-x+1/2, -y, z+1/2'
'-x, y+1/2, -z+1/2'
'x+1/2, -y+1/2, -z'
'-x, -y, -z'
'x-1/2, y, -z-1/2'
'x, -y-1/2, z-1/2'
'-x-1/2, y-1/2, z'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_calc_flag
_atom_site_occupancy
_atom_site_type_symbol
C1 -0.0537(8) 0.1425(9) 0.0097(12) 0.037(2) Uani d 1 C
H1 -0.085(6) 0.246(7) 0.034(8) 0.044 Uiso d 1 H
C2 0.0840(7) 0.0924(10) 0.1373(10) 0.040(2) Uani d 1 C
H2 0.140(6) 0.156(6) 0.219(8) 0.048 Uiso d 1 H
C3 0.1343(7) -0.0521(9) 0.1235(12) 0.044(2) Uani d 1 C
H3 0.220(6) -0.080(6) 0.204(9) 0.052 Uiso d 1 H
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
C1 0.039(3) 0.028(8) 0.043(6) 0.005(4) 0.002(4) 0.002(5)
C2 0.034(4) 0.060(7) 0.025(6) 0.005(4) -0.005(3) -0.001(4)
C3 0.043(4) 0.048(7) 0.039(6) 0.002(5) -0.007(4) 0.006(5)
loop_
_atom_type_symbol
_atom_type_description
_atom_type_scat_dispersion_real
_atom_type_scat_dispersion_imag
_atom_type_scat_source
C C 0.0033 0.0016 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
H H 0.0000 0.0000 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
loop_
_geom_angle_atom_site_label_1
_geom_angle_atom_site_label_2
_geom_angle_atom_site_label_3
_geom_angle_site_symmetry_1
_geom_angle
C3 C1 C2 5 120.8(8)
C3 C1 H1 5 127(3)
C2 C1 H1 . 113(3)
C3 C2 C1 . 117.9(7)
C3 C2 H2 . 123(4)
C1 C2 H2 . 119(4)
C1 C3 C2 5 121.4(7)
C1 C3 H3 5 123(4)
C2 C3 H3 . 116(4)
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_site_symmetry_2
_geom_bond_distance
C1 C3 5 1.353(10)
C1 C2 . 1.396(8)
C1 H1 . 0.99(6)
C2 C3 . 1.382(10)
C2 H2 . 0.90(5)
C3 C1 5 1.353(10)
C3 H3 . 0.86(5)
loop_
_geom_torsion_atom_site_label_1
_geom_torsion_atom_site_label_2
_geom_torsion_atom_site_label_3
_geom_torsion_atom_site_label_4
_geom_torsion_site_symmetry_1
_geom_torsion_site_symmetry_4
_geom_torsion
C3 C1 C2 C3 5 . -0.2(10)
C1 C2 C3 C1 . 5 0.2(10)
""")
    assign, lattice_info = Xponge.get_assignment_from_cif(cif)
    assign.determine_atom_type("gaff")
    assign.calculate_charge("tpacm4")
    restype = assign.to_residuetype("SOM")
    gaff.parmchk2_gaff(restype, 'som.frcmod')
    lattice = Xponge.Lattice(basis_molecule = restype, **lattice_info)
    region = Xponge.BlockRegion(0, 0, 0,
        5 * lattice_info["cell_length"][0],
        3 * lattice_info["cell_length"][1],
        3 * lattice_info["cell_length"][2])
    box = Xponge.BlockRegion(0, 0, 0,
        6 * lattice_info["cell_length"][0],
        4 * lattice_info["cell_length"][1],
        4 * lattice_info["cell_length"][2])
    mol = lattice.create(box, region)
    Xponge.save_mol2(mol, "som.mol2")
    Xponge.save_sponge_input(mol, "som")
    #assert run("SPONGE -default_in_file_prefix som -mode NPT -dt 1e-3 \
#-thermostat middle_langevin -barostat andersen_barostat -step_limit 5000 > som.out") == 0
