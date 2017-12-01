# gen_read_callsite_file.py

from parser import statements, block_statements, typedecl_statements
from kgplugin import Kgen_Plugin
from verify_utils import VERIFY_PBLOCK_USE_PART, VERIFY_PBLOCK_DECL_PART, VERIFY_PBLOCK_EXEC_PART, \
    VERIFY_PBLOCK_CONTAINS_PART, VERIFY_PBLOCK_SUBP_PART, VERIFY_PBLOCK_EXTERNS, VERIFY_PBLOCK_LOCALS, \
    VERIFY_PBLOCK_INIT

class Verify_K_Callsite_File(Kgen_Plugin):
    def __init__(self):
        self.frame_msg = None

    # registration
    def register(self, msg):
        self.frame_msg = msg

        # register initial events

        self.frame_msg.add_event(KERNEL_SELECTION.ALL, FILE_TYPE.KERNEL, GENERATION_STAGE.NODE_CREATED, \
            getinfo('parentblock_stmt'), None, self.create_parentblock_parts)

        self.frame_msg.add_event(KERNEL_SELECTION.ALL, FILE_TYPE.KERNEL, GENERATION_STAGE.NODE_CREATED, \
            getinfo('topblock_stmt'), None, self.create_topblock_parts)

    def create_topblock_parts(self, node):

        attrs = {'name': 'kgen_utils_mod', 'isonly': True, 'items':['check_t', 'kgen_init_check', 'kgen_tolerance', \
            'kgen_minvalue', 'CHECK_IDENTICAL', 'CHECK_IN_TOL', 'CHECK_OUT_TOL']}
        part_append_genknode(node, USE_PART, statements.Use, attrs=attrs)

        prenode = getinfo('blocknode_aftercallsite_main')
        self.frame_msg.add_event(KERNEL_SELECTION.ALL, FILE_TYPE.KERNEL, GENERATION_STAGE.BEGIN_PROCESS, \
            prenode, None, self.create_verification_parts)

        namedpart_create_subpart(prenode, VERIFY_PBLOCK_INIT, EXEC_PART)
        namedpart_append_comment(prenode.kgen_kernel_id, VERIFY_PBLOCK_INIT, '')
        namedpart_append_comment(prenode.kgen_kernel_id, VERIFY_PBLOCK_INIT, 'verify init')

        namedpart_create_subpart(prenode, VERIFY_PBLOCK_EXTERNS, EXEC_PART)
        namedpart_append_comment(prenode.kgen_kernel_id, VERIFY_PBLOCK_EXTERNS, '')
        namedpart_append_comment(prenode.kgen_kernel_id, VERIFY_PBLOCK_EXTERNS, 'extern verify variables')

        namedpart_create_subpart(prenode, VERIFY_PBLOCK_LOCALS, EXEC_PART)
        namedpart_append_comment(prenode.kgen_kernel_id, VERIFY_PBLOCK_LOCALS, '')
        namedpart_append_comment(prenode.kgen_kernel_id, VERIFY_PBLOCK_LOCALS, 'local verify variables')

    def create_verification_parts(self, node):

        attrs = {'designator': 'kgen_init_check', 'items': ['check_status', 'tolerance=%s'%getinfo('verify_tol'), \
            'verboseLevel=%s'%getinfo('verbose_level')]}
        namedpart_append_genknode(node.kgen_kernel_id, VERIFY_PBLOCK_INIT, statements.Call, attrs=attrs)

        attrs = {'items': ['""']}
        part_append_genknode(node, EXEC_PART, statements.Write, attrs=attrs)

        # verification statistics
        attrs = {'expr': 'check_status%verboseLevel > 0'}
        ifstatobj = part_append_genknode(node, EXEC_PART, block_statements.IfThen, attrs=attrs)

        attrs = {'items': ['"Number of output variables: "','check_status%numTotal']}
        part_append_genknode(ifstatobj, EXEC_PART, statements.Write, attrs=attrs)

        attrs = {'items': ['"Number of identical variables: "','check_status%numIdentical']}
        part_append_genknode(ifstatobj, EXEC_PART, statements.Write, attrs=attrs)

        attrs = {'items': ['"Number of non-identical variables within tolerance: "','check_status%numInTol']}
        part_append_genknode(ifstatobj, EXEC_PART, statements.Write, attrs=attrs)

        attrs = {'items': ['"Number of non-identical variables out of tolerance: "','check_status%numOutTol']}
        part_append_genknode(ifstatobj, EXEC_PART, statements.Write, attrs=attrs)

        attrs = {'items': ['"Tolerance: "','kgen_tolerance']}
        part_append_genknode(ifstatobj, EXEC_PART, statements.Write, attrs=attrs)

        attrs = {'items': ['""']}
        part_append_genknode(node, EXEC_PART, statements.Write, attrs=attrs)

        #jgw# do mpi_allreduce on numOutTol
        part_append_comment(node, EXEC_PART, '')
        attrs = {'variable': 'send(1)', 'sign': '=', 'expr': 'check_status%numOutTol'}
        part_append_genknode(node, EXEC_PART, statements.Assignment, attrs=attrs)
        attrs = {'designator': 'mpi_allreduce', 'items': ['send','recv','1','MPI_INT','MPI_MAX', getinfo('mpi_comm'), 'kgen_ierr']}
        part_append_gensnode(node, EXEC_PART, statements.Call, attrs=attrs)

        attrs = {'variable': 'check_status%numOutTol', 'sign': '=', 'expr': 'recv(1)'}
        part_append_genknode(node, EXEC_PART, statements.Assignment, attrs=attrs)

        attrs = {'designator': 'mpi_comm_rank', 'items': [getinfo('mpi_comm'),\
 'kgen_mpirank', 'kgen_ierr']}
        part_append_gensnode(node, EXEC_PART, statements.Call, attrs=attrs)
        part_append_comment(node, EXEC_PART, '')
        # jgw

        # verification result
        attrs = {'expr': 'check_status%numOutTol > 0'}
        ifobj = part_append_genknode(node, EXEC_PART, block_statements.IfThen, attrs=attrs)

        #jgw# only print verified on rank 0
        attrs = {'expr': 'kgen_mpirank == 0'}
        ifzero = part_append_genknode(ifobj, EXEC_PART, block_statements.IfThen, attrs=attrs)

        attrs = {'items': ['"Verification FAILED"']}
        #jgw#
        part_append_genknode(ifzero, EXEC_PART, statements.Write, attrs=attrs)
        
        attrs = {'variable': 'check_status%Passed', 'sign': '=', 'expr': '.FALSE.'}
        part_append_genknode(ifobj, EXEC_PART, statements.Assignment, attrs=attrs)

        attrs = {'variable': 'kgen_isverified', 'sign': '=', 'expr': '.FALSE.'}
        part_append_genknode(ifobj, EXEC_PART, statements.Assignment, attrs=attrs)

        part_append_genknode(ifobj, EXEC_PART, statements.Else, attrs=attrs)

        #jgw# only print verified on rank 0
        attrs = {'expr': 'kgen_mpirank == 0'}
        ifzero = part_append_genknode(ifobj, EXEC_PART, block_statements.IfThen, attrs=attrs)

        attrs = {'items': ['"Verification PASSED"']}
        #jgw#
        part_append_genknode(ifzero, EXEC_PART, statements.Write, attrs=attrs)
        
        attrs = {'variable': 'check_status%Passed', 'sign': '=', 'expr': '.TRUE.'}
        part_append_genknode(ifobj, EXEC_PART, statements.Assignment, attrs=attrs)

        attrs = {'variable': 'kgen_isverified', 'sign': '=', 'expr': '.TRUE.'}
        part_append_genknode(ifobj, EXEC_PART, statements.Assignment, attrs=attrs)

        attrs = {'items': ['""']}
        part_append_genknode(node, EXEC_PART, statements.Write, attrs=attrs)

    def create_parentblock_parts(self, node):

        namedpart_link_part(node, VERIFY_PBLOCK_USE_PART, USE_PART)
        namedpart_link_part(node, VERIFY_PBLOCK_DECL_PART, DECL_PART)
        namedpart_link_part(node, VERIFY_PBLOCK_EXEC_PART, EXEC_PART)
        namedpart_link_part(node, VERIFY_PBLOCK_CONTAINS_PART, CONTAINS_PART)
        namedpart_link_part(node, VERIFY_PBLOCK_SUBP_PART, SUBP_PART)

        #jgw# add use mpi. This was the best place I could find to put it,
        # but it doesn't seem quite right.
        attrs = {'name':'mpi', 'isonly': False, 'items':[]}
        part_append_genknode(node, USE_PART, statements.Use, attrs=attrs)

        attrs = {'type_spec': 'TYPE', 'selector':(None, 'check_t'), 'entity_decls': ['check_status']}
        part_append_genknode(node, DECL_PART, typedecl_statements.Type, attrs=attrs)
