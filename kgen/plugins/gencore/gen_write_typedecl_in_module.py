# gen_write_typedecl_in_module.py
 
import statements
import block_statements
import typedecl_statements
from kgen_plugin import Kgen_Plugin

from gencore_utils import get_topname, get_module_writename, MODULE_EXTERNS_PART, PARENTBLOCK_WRITE_IN_EXTERNS, \
    PARENTBLOCK_USE_PART

class Gen_S_Typedecl_In_Module(Kgen_Plugin):
    def __init__(self):
        self.frame_msg = None
        self.externs_part = None
        self.subrname = None

    # registration
    def register(self, msg):
        self.frame_msg = msg

        # register initial events
        self.frame_msg.add_event(KERNEL_SELECTION.ALL, FILE_TYPE.STATE, GENERATION_STAGE.NODE_CREATED, \
            block_statements.Module, self.has_externs_in_module, self.create_module_parts) 

    def has_externs_in_module(self, node):
        checks = lambda n: hasattr(n.kgen_stmt, 'geninfo') and len(n.kgen_stmt.geninfo)>0 \
            and isinstance(n.kgen_stmt, typedecl_statements.TypeDeclarationStatement) \
            and "parameter" not in n.kgen_stmt.attrspec
        if part_has_node(node, DECL_PART, checks):
            return True
        return False

    def is_extern_in_module(self, node):
        if node.kgen_stmt and hasattr(node.kgen_stmt, 'geninfo') and len(node.kgen_stmt.geninfo)>0 and \
            'parameter' not in node.kgen_stmt.attrspec:
            return True
        return False

    def create_module_parts(self, node):
        node.kgen_stmt.top.used4genstate = True

        self.modname = node.name
        self.subrname = get_module_writename(node.kgen_stmt)
        attrs = {'name': self.subrname, 'args': ['kgen_unit']}
        subrobj = part_append_gensnode(node, SUBP_PART, block_statements.Subroutine, attrs=attrs)

        # kgen_unit
        attrs = {'type_spec': 'INTEGER', 'attrspec': ['INTENT(IN)'], 'entity_decls': ['kgen_unit']}
        part_append_gensnode(subrobj, DECL_PART, typedecl_statements.Integer, attrs=attrs)

        # add public stmt
        attrs = {'items':[self.subrname]}
        part_append_gensnode(node, DECL_PART, statements.Public, attrs=attrs)

        namedpart_link_part(subrobj, MODULE_EXTERNS_PART, EXEC_PART)

        # register event per typedecl 
        self.frame_msg.add_event(KERNEL_SELECTION.ALL, FILE_TYPE.STATE, GENERATION_STAGE.BEGIN_PROCESS, \
            typedecl_statements.TypeDeclarationStatement, self.is_extern_in_module, self.create_subr_write_typedecl_in_module) 

        # register event per module
        self.frame_msg.add_event(KERNEL_SELECTION.ALL, FILE_TYPE.STATE, GENERATION_STAGE.BEGIN_PROCESS, \
            block_statements.Module, self.has_externs_in_module, self.create_stmts_in_callsite) 


    def create_stmts_in_callsite(self, node):
        checks = lambda n: isinstance(n.kgen_stmt, statements.Use) and n.kgen_stmt.name==self.modname and \
            ( not n.kgen_stmt.isonly or self.subrname in [ item.split('=>')[0].strip() for item in n.kgen_stmt.items])
        if not namedpart_has_node(node.kgen_kernel_id, PARENTBLOCK_USE_PART, checks):
            attrs = {'name':self.modname, 'isonly': True, 'items':[self.subrname]}
            namedpart_append_gensnode(node.kgen_kernel_id, PARENTBLOCK_USE_PART, statements.Use, attrs=attrs)

        attrs = {'designator': self.subrname, 'items': ['kgen_unit']}
        namedpart_append_gensnode(node.kgen_kernel_id, PARENTBLOCK_WRITE_IN_EXTERNS, statements.Call, attrs=attrs)

    def create_subr_write_typedecl_in_module(self, node):
        pass
        #import pdb; pdb.set_trace()
        # create write or call stmt

        # create subr if necessary


#
#    def create_write_subr(self, subrname, entity_name, parent, var, stmt):
#        if subrname not in self.created_subrs and \
#            not has_object_in_part(parent, SUBP_PART, block_statements.Subroutine, attrs={'name':subrname} ):
#
#            self.created_subrs.append(subrname)
#
#            if not self.is_contains_created and not has_object_in_part(parent, CONTAINS_PART, statements.Contains ):
#                append_comment_in_part(parent, CONTAINS_PART, '')
#                part_append_gensnode(parent, CONTAINS_PART, gensobj(parent, statements.Contains, parent.kgen_kernel_id))
#                append_comment_in_part(parent, CONTAINS_PART, '')
#                self.is_contains_created = True
#
#            append_comment_in_part(parent, SUBP_PART, 'write state subroutine for %s'%subrname)
#            attrs = {'name': subrname, 'args': ['var', 'kgen_unit', 'printvar']}
#            subrobj = gensobj(parent, block_statements.Subroutine, parent.kgen_kernel_id, attrs=attrs)
#            part_append_gensnode(parent, SUBP_PART, subrobj)
#            append_comment_in_part(parent, SUBP_PART, '')
#
#            parent.kgen_stmt.top.used4genstate = True
#
#            # variable A
#            #import pdb; pdb.set_trace()
#            attrspec = ['INTENT(IN)']
#            #if var.is_pointer(): attrspec.append('POINTER')
#            if var.is_array(): attrspec.append('DIMENSION(%s)'% ','.join(':'*var.rank))
#            attrs = {'type_spec': stmt.name.upper(), 'attrspec': attrspec, 'selector':stmt.selector, 'entity_decls': ['var']}
#            part_append_gensnode(subrobj, DECL_PART, gensobj(subrobj, stmt.__class__, subrobj.kgen_kernel_id, attrs=attrs))
#
#            # kgen_unit
#            attrs = {'type_spec': 'INTEGER', 'attrspec': ['INTENT(IN)'], 'entity_decls': ['kgen_unit']}
#            part_append_gensnode(subrobj, DECL_PART, gensobj(subrobj, typedecl_statements.Integer, subrobj.kgen_kernel_id, attrs=attrs))
#
#            # printvar
#            attrs = {'type_spec': 'CHARACTER', 'attrspec': ['INTENT(IN)', 'OPTIONAL'], 'selector':('*', None), 'entity_decls': ['printvar']}
#            part_append_gensnode(subrobj, DECL_PART, gensobj(subrobj, typedecl_statements.Character, subrobj.kgen_kernel_id, attrs=attrs))
#
#            # is_true
#            attrs = {'type_spec': 'LOGICAL', 'entity_decls': ['is_true']}
#            part_append_gensnode(subrobj, DECL_PART, gensobj(subrobj, typedecl_statements.Logical, subrobj.kgen_kernel_id, attrs=attrs))
#
#            # array index A
#            if var.is_array():
#                attrs = {'type_spec': 'INTEGER', 'entity_decls': [ 'idx%d'%(d+1) for d in range(var.rank) ]}
#                part_append_gensnode(subrobj, DECL_PART, gensobj(subrobj, typedecl_statements.Integer, subrobj.kgen_kernel_id, attrs=attrs))
#
#            append_comment_in_part(subrobj, DECL_PART, '')
#
#            #ifobj = None
#
##            # if var is pointer
##            if var.is_pointer():
##                attrs = {'expr': '.NOT. ASSOCIATED(var)'}
##                ifptrobj = gensobj(subrobj, block_statements.IfThen, subrobj.kgen_kernel_id, attrs=attrs)
##                part_append_gensnode(subrobj, EXEC_PART, ifptrobj)
##
##                ifobj = ifptrobj
##
##                attrs = {'variable': 'is_true', 'sign': '=', 'expr': '.FALSE.'}
##                part_append_gensnode(ifptrobj, EXEC_PART, gensobj(ifptrobj, statements.Assignment, ifptrobj.kgen_kernel_id, attrs=attrs))
##
##                # if var is array
##                if var.is_array():
##                    attrs = {'expr': 'SIZE(var)==1'}
##                    part_append_gensnode(ifptrobj, EXEC_PART, gensobj(ifptrobj, block_statements.ElseIf, ifptrobj.kgen_kernel_id, attrs=attrs))
##
##            # if var is array
##            elif var.is_array():
#            attrs = {'expr': 'SIZE(var)==1'}
#            ifarrobj = gensobj(subrobj, block_statements.IfThen, subrobj.kgen_kernel_id, attrs=attrs)
#            part_append_gensnode(subrobj, EXEC_PART, ifarrobj)
#
#            ifobj = ifarrobj
#
#            if var.is_array():
#                attrs = {'expr': 'UBOUND(var, 1)<LBOUND(var, 1)'}
#                ifbndobj = gensobj(ifobj, block_statements.IfThen, ifobj.kgen_kernel_id, attrs=attrs)
#                part_append_gensnode(ifobj, EXEC_PART, ifbndobj)
#
#                attrs = {'variable': 'is_true', 'sign': '=', 'expr': '.FALSE.'}
#                part_append_gensnode(ifbndobj, EXEC_PART, gensobj(ifbndobj, statements.Assignment, ifbndobj.kgen_kernel_id, attrs=attrs))
#
#                attrs = {'expr': 'UBOUND(var, 1)==0 .AND. LBOUND(var, 1)==0'}
#                part_append_gensnode(ifbndobj, EXEC_PART, gensobj(ifbndobj, block_statements.ElseIf, ifbndobj.kgen_kernel_id, attrs=attrs))
#
#                attrs = {'variable': 'is_true', 'sign': '=', 'expr': '.FALSE.'}
#                part_append_gensnode(ifbndobj, EXEC_PART, gensobj(ifbndobj, statements.Assignment, ifbndobj.kgen_kernel_id, attrs=attrs))
#
#                part_append_gensnode(ifbndobj, EXEC_PART, gensobj(ifbndobj, block_statements.Else, ifbndobj.kgen_kernel_id, attrs=attrs))
#
#                attrs = {'variable': 'is_true', 'sign': '=', 'expr': '.TRUE.'}
#                part_append_gensnode(ifbndobj, EXEC_PART, gensobj(ifbndobj, statements.Assignment, ifbndobj.kgen_kernel_id, attrs=attrs))
#
#            part_append_gensnode(ifobj, EXEC_PART, gensobj(ifobj, block_statements.Else, ifobj.kgen_kernel_id, attrs=attrs))
#
#            attrs = {'variable': 'is_true', 'sign': '=', 'expr': '.TRUE.'}
#            part_append_gensnode(ifobj, EXEC_PART, gensobj(ifobj, statements.Assignment, ifobj.kgen_kernel_id, attrs=attrs))
#            append_comment_in_part(subrobj, EXEC_PART, '')
#
#            attrs = {'items': ['is_true'], 'specs': ['UNIT = kgen_unit']}
#            part_append_gensnode(subrobj, EXEC_PART, gensobj(subrobj, statements.Write, subrobj.kgen_kernel_id, attrs=attrs))
#            append_comment_in_part(subrobj, EXEC_PART, '')
#
#            attrs = {'expr': 'is_true'}
#            iftrueobj = gensobj(subrobj, block_statements.IfThen, subrobj.kgen_kernel_id, attrs=attrs)
#            part_append_gensnode(subrobj, EXEC_PART, iftrueobj)
#
#            if var.is_array():
#                for dim in range(var.rank):
#                    attrs = {'items': ['LBOUND(var, %d)'%(dim+1)], 'specs': ['UNIT = kgen_unit']}
#                    part_append_gensnode(iftrueobj, EXEC_PART, gensobj(iftrueobj, statements.Write, iftrueobj.kgen_kernel_id, attrs=attrs))
#
#                    attrs = {'items': ['UBOUND(var, %d)'%(dim+1)], 'specs': ['UNIT = kgen_unit']}
#                    part_append_gensnode(iftrueobj, EXEC_PART, gensobj(iftrueobj, statements.Write, iftrueobj.kgen_kernel_id, attrs=attrs))
#
#                if stmt.is_derived():
#                    indexes = [ 'idx%d'%(d+1) for d in range(var.rank) ]
#                    str_indexes = ','.join(indexes)
#
#                    prevobj = iftrueobj
#                    doobjs = []
#                    for d in range(var.rank):
#                        attrs = {'loopcontrol': 'idx%(d)d=LBOUND(var,%(d)d), UBOUND(var,%(d)d)'%{'d':d+1}}
#                        doobj = gensobj(prevobj, block_statements.Do, prevobj.kgen_kernel_id, attrs=attrs)
#                        part_append_gensnode(prevobj, EXEC_PART, doobj)
#
#                        doobjs.append(doobj)
#                        prevobj = doobj
#
#                    attrs = {'expr': 'PRESENT( printvar )'}
#                    ifpvarobj = gensobj(doobjs[-1], block_statements.IfThen, doobjs[-1].kgen_kernel_id, attrs=attrs)
#                    part_append_gensnode(doobjs[-1], EXEC_PART, ifpvarobj)
#
#                    callname = None
#                    for uname, req in stmt.unknowns.iteritems():
#                        if uname.firstpartname()==stmt.name:
#                            res = req.res_stmts[0]
#                            callname = get_dtype_writename(res)
#                            break
#                    if callname is None: raise ProgramException('Can not find Type resolver for %s'%stmt.name)
#
#                    attrs = {'designator': callname, 'items': ['var(%s)'%str_indexes, 'kgen_unit', 'printvar // "(%s)"'%str_indexes]}
#                    part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Call, ifpvarobj.kgen_kernel_id, attrs=attrs))
#
#                    part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Else, ifpvarobj.kgen_kernel_id))
#
#                    if any(match_namepath(pattern, pack_exnamepath(stmt, entity_name), internal=False) for pattern in getinfo('print_var_names')):
#                        attrs = {'designator': callname, 'items': ['var(%s)'%str_indexes, 'kgen_unit', '"%s(%s)"'%(entity_name, str_indexes)]}
#                        part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Call, ifpvarobj.kgen_kernel_id, attrs=attrs))
#                    else:
#                        attrs = {'designator': callname, 'items': ['var%%%s'%entity_name, 'kgen_unit']}
#                        part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Call, ifpvarobj.kgen_kernel_id, attrs=attrs))
#
#                else: # intrinsic type
#                    attrs = {'items': ['var'], 'specs': ['UNIT = kgen_unit']}
#                    part_append_gensnode(iftrueobj, EXEC_PART, gensobj(iftrueobj, statements.Write, iftrueobj.kgen_kernel_id, attrs=attrs))
#
#                    if any(match_namepath(pattern, pack_exnamepath(stmt, entity_name), internal=False) for pattern in getinfo('print_var_names')):
#                        attrs = {'items': ['"** KGEN DEBUG: " // "%s **"'%entity_name, 'var']}
#                        part_append_gensnode(iftrueobj, EXEC_PART, gensobj(iftrueobj, statements.Write, iftrueobj.kgen_kernel_id, attrs=attrs))
#                    else:
#                        attrs = {'expr': 'PRESENT( printvar )'}
#                        ifpvarobj = gensobj(iftrueobj, block_statements.IfThen, iftrueobj.kgen_kernel_id, attrs=attrs)
#                        part_append_gensnode(iftrueobj, EXEC_PART, ifpvarobj)
#
#                        attrs = {'items': ['"** KGEN DEBUG: " // printvar // " %s **"'%entity_name, 'var']}
#                        part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Write, ifpvarobj.kgen_kernel_id, attrs=attrs))
#
#
#            else: # scalar
#                if self.stmt.is_derived():
#                    attrs = {'expr': 'PRESENT( printvar )'}
#                    ifpvarobj = gensobj(iftrueobj, block_statements.IfThen, iftrueobj.kgen_kernel_id, attrs=attrs)
#                    part_append_gensnode(iftrueobj, EXEC_PART, ifpvarobj)
#
#                    callname = None
#                    for uname, req in stmt.unknowns.iteritems():
#                        if uname.firstpartname()==stmt.name:
#                            res = req.res_stmts[0]
#                            callname = get_dtype_writename(res)
#                            break
#                    if callname is None: raise ProgramException('Can not find Type resolver for %s'%stmt.name)
#
#                    attrs = {'designator': callname, 'items': ['var', 'kgen_unit', 'printvar // " %s "'%entity_name]}
#                    part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Call, ifpvarobj.kgen_kernel_id, attrs=attrs))
#
#                    part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Else, ifpvarobj.kgen_kernel_id))
#
#                    if any(match_namepath(pattern, pack_exnamepath(stmt, entity_name), internal=False) for pattern in getinfo('print_var_names')):
#                        attrs = {'designator': callname, 'items': ['var', 'kgen_unit', '"%s"'%entity_name]}
#                        part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Call, ifpvarobj.kgen_kernel_id, attrs=attrs))
#                    else:
#                        attrs = {'designator': callname, 'items': ['var', 'kgen_unit']}
#                        part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Call, ifpvarobj.kgen_kernel_id, attrs=attrs))
#                else: # intrinsic type
#                    attrs = {'items': ['var'], 'specs': ['UNIT = kgen_unit']}
#                    part_append_gensnode(iftrueobj, EXEC_PART, gensobj(iftrueobj, statements.Write, iftrueobj.kgen_kernel_id, attrs=attrs))
#
#                    if any(match_namepath(pattern, pack_exnamepath(stmt, entity_name), internal=False) for pattern in getinfo('print_var_names')):
#                        attrs = {'items': ['"** KGEN DEBUG: " // "%s **"'%entity_name, 'var']}
#                        part_append_gensnode(iftrueobj, EXEC_PART, gensobj(iftrueobj, statements.Write, iftrueobj.kgen_kernel_id, attrs=attrs))
#                    else:
#                        attrs = {'expr': 'PRESENT( printvar )'}
#                        ifpvarobj = gensobj(iftrueobj, block_statements.IfThen, iftrueobj.kgen_kernel_id, attrs=attrs)
#                        part_append_gensnode(iftrueobj, EXEC_PART, ifpvarobj)
#
#                        attrs = {'items': ['"** KGEN DEBUG: " // printvar // " %s **"'%entity_name, 'var']}
#                        part_append_gensnode(ifpvarobj, EXEC_PART, gensobj(ifpvarobj, statements.Write, ifpvarobj.kgen_kernel_id, attrs=attrs))
#
#       
#    def create_subr_write_typedecl_in_type(self, node):
#
#        stmt = node.kgen_stmt
#        parent = node.kgen_parent.kgen_parent
#        entity_names = [ get_entity_name(decl) for decl in stmt.entity_decls ]
#
#        for entity_name, entity_decl in zip(entity_names, stmt.entity_decls):
#            var = stmt.get_variable(entity_name)
#            subrname = get_typedecl_subpname(stmt, entity_name)
#            if subrname is None: raise ProgramException('Can not get subroutinename')
#
#            if var.is_array():
#                if stmt.is_derived():
#                    self.create_write_subr(subrname, entity_name, parent, var, stmt)
#                else: # intrinsic type
#                    if var.is_explicit_shape_array():
#                        pass
#                    else: # implicit array
#                        self.create_write_subr(subrname, entity_name, parent, var, stmt)
#
