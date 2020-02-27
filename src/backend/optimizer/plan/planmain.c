/*-------------------------------------------------------------------------
 *
 * planmain.c
 *	  Routines to plan a single query
 *
 * What's in a name, anyway?  The top-level entry point of the planner/
 * optimizer is over in planner.c, not here as you might think from the
 * file name.  But this is the main code for planning a basic join operation,
 * shorn of features like subselects, inheritance, aggregates, grouping,
 * and so on.  (Those are the things planner.c deals with.)
 *
 * Portions Copyright (c) 1996-2019, PostgreSQL Global Development Group
 * Portions Copyright (c) 1994, Regents of the University of California
 *
 *
 * IDENTIFICATION
 *	  src/backend/optimizer/plan/planmain.c
 *
 *-------------------------------------------------------------------------
 */
#include "postgres.h"

#include "optimizer/appendinfo.h"
#include "optimizer/clauses.h"
#include "optimizer/inherit.h"
#include "optimizer/optimizer.h"
#include "optimizer/orclauses.h"
#include "optimizer/pathnode.h"
#include "optimizer/paths.h"
#include "optimizer/placeholder.h"
#include "optimizer/planmain.h"

int
generate_valid_join_orders_recursive(PlannerInfo * root, struct JoinTreeNode * node, List ** children, Bitmapset ** nodeBitmapSets, List ** validBitmapSets,
	int * visitOrder, int * count)
{
	int level = 0;
	Bitmapset * bitmap = NULL;
	printf("[generate_valid_join_orders_recursive] Start processing node %d\n", node->nodeId);
	// Process base relation.
	if (children[node->nodeId] == NULL) {
		ListCell * lc2;

		// Check the corresponding relid for this node.
		bool flag = false;
		printf("[generate_valid_join_orders_recursive] Number of range table entries: %d\n", root->simple_rel_array_size - 1);
		for (int rti = 1; rti < root->simple_rel_array_size; ++rti) {
			RangeTblEntry * rte = root->simple_rte_array[rti];
			if (rte == NULL) {
				continue;
			}
			printf("[generate_valid_join_orders_recursive] Range table entry: [name %s, alias %s, relid %d, index %d]\n",
				rte->eref->aliasname, rte->alias == NULL ? "null" : rte->alias->aliasname, rte->relid,
					rti);
			if (strcmp(rte->eref->aliasname, node->relName) == 0
				&& (node->alias == NULL || strcmp(rte->alias->aliasname, node->alias) == 0))
			{
				printf("[generate_valid_join_orders_recursive] Found range table entry [name %s, alias %s, relid %d, index %d] for node [id %d, name %s, alias %s], visit order %d\n",
					rte->eref->aliasname, rte->alias == NULL ? "null" : rte->alias->aliasname, rte->relid, rti,
					node->nodeId, node->relName, node->alias == NULL ? "null" : node->alias, *count);
				visitOrder[rti] = *count;
				*count += 1;
				bitmap = bms_add_member(bitmap, rti);
				flag = true;
				break;
			}
		}
		// Do not force the plan if a table is not found in the range table entries.
		if (!flag)
		{
			printf("[generate_valid_join_orders_recursive] Cannot force the plan because relation [%s, %s] is not found in range tables\n",
				node->relName, node->alias == NULL ? "null" : node->alias);
			return -1;
		}
		level = 1;
	}
	else
	{
		ListCell * lc;
		foreach(lc, children[node->nodeId]) {
			struct JoinTreeNode * child = (struct JoinTreeNode *)lfirst(lc);
			int m = generate_valid_join_orders_recursive(root, child, children, nodeBitmapSets, validBitmapSets,
				visitOrder, count);
			if (m == -1) {
				// There exists a base relation in the subtree that cannot be found.
				return -1;
			}
			level += m;
			printf("[generate_valid_join_orders_recursive] bitmap before %s\n",
					bitmap == NULL ? "null" : bmsToString(bitmap));
			bitmap = bms_add_members(bitmap, nodeBitmapSets[child->nodeId]);
			printf("[generate_valid_join_orders_recursive] Add bitmap, child %d, parent %d\n",
				child->nodeId, node->nodeId);
		}
	}
	nodeBitmapSets[node->nodeId] = bitmap;
	validBitmapSets[level] = lappend(validBitmapSets[level], bitmap);
	printf("[Bailu.DEBUG]: Finish processing node %d at level %d with bitmap %s\n",
		node->nodeId, level, bmsToString(bitmap));
	return level;
}

 /*
 * Bailu:
 * generate_valid_join_orders
 *  If join order is forced, generate all subsets of base relations that are compatible with
 *  the forced join order.
 */
void
generate_valid_join_orders(PlannerInfo * root, List * joinTreeNodeList)
{
	int n = list_length(joinTreeNodeList);

	// Generate children information.
	List ** children = (List **)palloc0(n * sizeof(List *));
	for (int i = 0; i < n; ++i) {
		children[i] = NULL;
	}
	ListCell * lc;
	struct JoinTreeNode * rootNode;
	foreach(lc, joinTreeNodeList) {
		struct JoinTreeNode * node = (struct JoinTreeNode *)lfirst(lc);
		struct JoinTreeNode * parent = node->parentNode;
		if (parent != NULL) {
			children[parent->nodeId] = lappend(children[parent->nodeId], node);
			printf("[generate_valid_join_orders] Add child %d to parent %d\n",
				node->nodeId, parent->nodeId);
		}
		else
		{
			rootNode = node;
		}
	}

	// Recursively generate the valid join orders.
	Bitmapset ** nodeBitmapSets = (Bitmapset **)palloc0(n * sizeof(Bitmapset *));
	List ** validBitmapSets = (List **)palloc0(n * sizeof(List *));
	int * visitOrder = (int *)palloc0(root->simple_rel_array_size * sizeof(int));
	int count = 0;
	for (int i = 0; i < n; ++i) {
		nodeBitmapSets[i] = NULL;
		validBitmapSets[i] = NULL;
	}
	int valid = generate_valid_join_orders_recursive(root, rootNode, children, nodeBitmapSets, validBitmapSets, visitOrder, &count);
	if (valid == -1)
	{
		root->glob->validJoinBitmapSets = NULL;
		root->glob->visitOrder = NULL;
		list_free_deep(validBitmapSets);
		pfree(visitOrder);
		return;
	}
	// Print stats of valid bitmap sets.
	for (int i = 1; i < n; ++i) {
		printf("[generate_valid_join_orders] # of valid join orders at level %d: %d\n",
			i, validBitmapSets[i] == NULL ? 0 : list_length(validBitmapSets[i]));
		ListCell * lc;
		foreach(lc, validBitmapSets[i]) {
			Bitmapset * bitmap = (Bitmapset *)lfirst(lc);
			printf("%s, ", bmsToString(bitmap));
		}
		printf("\n");
	}
	printf("[Bailu.DEBUG]: end of valid join orders\n");

	root->glob->validJoinBitmapSets = validBitmapSets;
	root->glob->visitOrder = visitOrder;
}

/*
 * query_planner
 *	  Generate a path (that is, a simplified plan) for a basic query,
 *	  which may involve joins but not any fancier features.
 *
 * Since query_planner does not handle the toplevel processing (grouping,
 * sorting, etc) it cannot select the best path by itself.  Instead, it
 * returns the RelOptInfo for the top level of joining, and the caller
 * (grouping_planner) can choose among the surviving paths for the rel.
 *
 * root describes the query to plan
 * qp_callback is a function to compute query_pathkeys once it's safe to do so
 * qp_extra is optional extra data to pass to qp_callback
 *
 * Note: the PlannerInfo node also includes a query_pathkeys field, which
 * tells query_planner the sort order that is desired in the final output
 * plan.  This value is *not* available at call time, but is computed by
 * qp_callback once we have completed merging the query's equivalence classes.
 * (We cannot construct canonical pathkeys until that's done.)
 */
RelOptInfo *
query_planner(PlannerInfo *root,
			  query_pathkeys_callback qp_callback, void *qp_extra)
{
	Query	   *parse = root->parse;
	List	   *joinlist;
	RelOptInfo *final_rel;

	/*
	 * Init planner lists to empty.
	 *
	 * NOTE: append_rel_list was set up by subquery_planner, so do not touch
	 * here.
	 */
	root->join_rel_list = NIL;
	root->join_rel_hash = NULL;
	root->join_rel_level = NULL;
	root->join_cur_level = 0;
	root->canon_pathkeys = NIL;
	root->left_join_clauses = NIL;
	root->right_join_clauses = NIL;
	root->full_join_clauses = NIL;
	root->join_info_list = NIL;
	root->placeholder_list = NIL;
	root->fkey_list = NIL;
	root->initial_rels = NIL;

	/*
	 * Make a flattened version of the rangetable for faster access (this is
	 * OK because the rangetable won't change any more), and set up an empty
	 * array for indexing base relations.
	 */
	setup_simple_rel_arrays(root);

	/*
	* Bailu
	* Process plan forcing after the rangetable has been set up.
	*/
	if (root->glob->joinTreeNodeList != NIL)
	{
		generate_valid_join_orders(root, root->glob->joinTreeNodeList);
		if (root->glob->validJoinBitmapSets == NIL) {
			// If no valid join order is constructed, disable plan forcing.
			root->glob->joinTreeNodeList = NIL;
		}
	}

	/*
	 * In the trivial case where the jointree is a single RTE_RESULT relation,
	 * bypass all the rest of this function and just make a RelOptInfo and its
	 * one access path.  This is worth optimizing because it applies for
	 * common cases like "SELECT expression" and "INSERT ... VALUES()".
	 */
	Assert(parse->jointree->fromlist != NIL);
	if (list_length(parse->jointree->fromlist) == 1)
	{
		Node	   *jtnode = (Node *) linitial(parse->jointree->fromlist);

		if (IsA(jtnode, RangeTblRef))
		{
			int			varno = ((RangeTblRef *) jtnode)->rtindex;
			RangeTblEntry *rte = root->simple_rte_array[varno];

			Assert(rte != NULL);
			if (rte->rtekind == RTE_RESULT)
			{
				/* Make the RelOptInfo for it directly */
				final_rel = build_simple_rel(root, varno, NULL);

				/*
				 * If query allows parallelism in general, check whether the
				 * quals are parallel-restricted.  (We need not check
				 * final_rel->reltarget because it's empty at this point.
				 * Anything parallel-restricted in the query tlist will be
				 * dealt with later.)  This is normally pretty silly, because
				 * a Result-only plan would never be interesting to
				 * parallelize.  However, if force_parallel_mode is on, then
				 * we want to execute the Result in a parallel worker if
				 * possible, so we must do this.
				 */
				if (root->glob->parallelModeOK &&
					force_parallel_mode != FORCE_PARALLEL_OFF)
					final_rel->consider_parallel =
						is_parallel_safe(root, parse->jointree->quals);

				/*
				 * The only path for it is a trivial Result path.  We cheat a
				 * bit here by using a GroupResultPath, because that way we
				 * can just jam the quals into it without preprocessing them.
				 * (But, if you hold your head at the right angle, a FROM-less
				 * SELECT is a kind of degenerate-grouping case, so it's not
				 * that much of a cheat.)
				 */
				add_path(final_rel, (Path *)
						 create_group_result_path(root, final_rel,
												  final_rel->reltarget,
												  (List *) parse->jointree->quals));

				/* Select cheapest path (pretty easy in this case...) */
				set_cheapest(final_rel);

				/*
				 * We still are required to call qp_callback, in case it's
				 * something like "SELECT 2+2 ORDER BY 1".
				 */
				(*qp_callback) (root, qp_extra);

				return final_rel;
			}
		}
	}

	/*
	 * Populate append_rel_array with each AppendRelInfo to allow direct
	 * lookups by child relid.
	 */
	setup_append_rel_array(root);

	/*
	 * Construct RelOptInfo nodes for all base relations used in the query.
	 * Appendrel member relations ("other rels") will be added later.
	 *
	 * Note: the reason we find the baserels by searching the jointree, rather
	 * than scanning the rangetable, is that the rangetable may contain RTEs
	 * for rels not actively part of the query, for example views.  We don't
	 * want to make RelOptInfos for them.
	 */
	add_base_rels_to_query(root, (Node *) parse->jointree);

	/*
	 * Examine the targetlist and join tree, adding entries to baserel
	 * targetlists for all referenced Vars, and generating PlaceHolderInfo
	 * entries for all referenced PlaceHolderVars.  Restrict and join clauses
	 * are added to appropriate lists belonging to the mentioned relations. We
	 * also build EquivalenceClasses for provably equivalent expressions. The
	 * SpecialJoinInfo list is also built to hold information about join order
	 * restrictions.  Finally, we form a target joinlist for make_one_rel() to
	 * work from.
	 */
	build_base_rel_tlists(root, root->processed_tlist);

	find_placeholders_in_jointree(root);

	find_lateral_references(root);

	joinlist = deconstruct_jointree(root);

	/*
	 * Reconsider any postponed outer-join quals now that we have built up
	 * equivalence classes.  (This could result in further additions or
	 * mergings of classes.)
	 */
	reconsider_outer_join_clauses(root);

	/*
	 * If we formed any equivalence classes, generate additional restriction
	 * clauses as appropriate.  (Implied join clauses are formed on-the-fly
	 * later.)
	 */
	generate_base_implied_equalities(root);

	/*
	 * We have completed merging equivalence sets, so it's now possible to
	 * generate pathkeys in canonical form; so compute query_pathkeys and
	 * other pathkeys fields in PlannerInfo.
	 */
	(*qp_callback) (root, qp_extra);

	/*
	 * Examine any "placeholder" expressions generated during subquery pullup.
	 * Make sure that the Vars they need are marked as needed at the relevant
	 * join level.  This must be done before join removal because it might
	 * cause Vars or placeholders to be needed above a join when they weren't
	 * so marked before.
	 */
	fix_placeholder_input_needed_levels(root);

	/*
	 * Remove any useless outer joins.  Ideally this would be done during
	 * jointree preprocessing, but the necessary information isn't available
	 * until we've built baserel data structures and classified qual clauses.
	 */
	joinlist = remove_useless_joins(root, joinlist);

	/*
	 * Also, reduce any semijoins with unique inner rels to plain inner joins.
	 * Likewise, this can't be done until now for lack of needed info.
	 */
	reduce_unique_semijoins(root);

	/*
	 * Now distribute "placeholders" to base rels as needed.  This has to be
	 * done after join removal because removal could change whether a
	 * placeholder is evaluable at a base rel.
	 */
	add_placeholders_to_base_rels(root);

	/*
	 * Construct the lateral reference sets now that we have finalized
	 * PlaceHolderVar eval levels.
	 */
	create_lateral_join_info(root);

	/*
	 * Match foreign keys to equivalence classes and join quals.  This must be
	 * done after finalizing equivalence classes, and it's useful to wait till
	 * after join removal so that we can skip processing foreign keys
	 * involving removed relations.
	 */
	match_foreign_keys_to_quals(root);

	/*
	 * Look for join OR clauses that we can extract single-relation
	 * restriction OR clauses from.
	 */
	extract_restriction_or_clauses(root);

	/*
	 * Now expand appendrels by adding "otherrels" for their children.  We
	 * delay this to the end so that we have as much information as possible
	 * available for each baserel, including all restriction clauses.  That
	 * let us prune away partitions that don't satisfy a restriction clause.
	 * Also note that some information such as lateral_relids is propagated
	 * from baserels to otherrels here, so we must have computed it already.
	 */
	add_other_rels_to_query(root);

	/*
	 * Ready to do the primary planning.
	 */
	final_rel = make_one_rel(root, joinlist);

	/* Check that we got at least one usable path */
	if (!final_rel || !final_rel->cheapest_total_path ||
		final_rel->cheapest_total_path->param_info != NULL)
		elog(ERROR, "failed to construct the join relation");

	return final_rel;
}
