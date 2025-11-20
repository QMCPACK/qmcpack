"""
EnhancedProjectManager class extending ProjectManager with enhanced DAG workflow support.

Provides conditional filtering, enhanced progress tracking, and blocked simulation
handling for DAG workflows with loops, conditionals, and error handlers.
"""

from typing import Dict, Set, List, Optional, Any
from developer import obj
from project_manager import ProjectManager, trivial
from simulation import Simulation
from enhanced_simulation import EnhancedSimulation


class EnhancedProjectManager(ProjectManager):
    """
    Extended ProjectManager with enhanced DAG workflow support.
    
    Inheritance Pattern:
    - Inherits from ProjectManager using standard Python inheritance
    - Calls super().__init__() to initialize base attributes (simulations, cascades, progressing_cascades)
    - Adds enhanced tracking: filtered_sims, blocked simulation reporting
    
    Method Overrides:
    
    1. add_simulations():
       Base: Adds all non-fake simulations to cascades
       Override: Adds pre-filtering for machine compatibility, definition-time conditionals,
                 and transitive dependency blocking. Snapshot blocked sims before detaching.
                 Calls super().add_simulations(final_verified) at end.
    
    2. run_project():
       Base: Initializes cascades and runs project
       Override: Detaches/prunes blocked sims BEFORE calling super().run_project()
    
    3. init_cascades():
       Base: Screens fake sims, resolves collisions, propagates blockages
       Override: Detaches/prunes blocked sims BEFORE calling super().init_cascades()
    
    4. screen_fake_sims():
       Base: Errors if any fake sims found
       Override: Allows intentionally blocked EnhancedSimulations (skips them in error check)
    
    5. traverse_cascades():
       Base: Standard DAG traversal
       Override: Uses standard DAG traversal (all workflows are DAGs)
    
    6. write_simulation_status():
       Base: Writes status for all simulations
       Override: Filters blocked sims from main status, adds separate "blocked simulations" section
    
    7. status_line():
       Base: Writes single status line
       Override: Adds [iter=N] annotation for simulations with loops (iteration_count > 0)
    
    New Methods (Not in Base):
    - _detach_blocked_simulations(): Removes blocked sims from dependency graph
    - _prune_blocked_from_tracking(): Cleans up tracking dictionaries
    - _snapshot_blocked(): Captures blocked sim info before graph mutation
    - _get_blocked_simulations(): Identifies all blocked sims (direct + transitive)
    
    Extension Strategy:
    - Pre-processing: Filters/validates before calling parent methods
    - Post-processing: Adds tracking/reporting after parent methods
    - Side-by-side: Adds new functionality (blocked sim reporting, enhanced status)
    - Selective override: Only overrides methods needed for enhanced DAG workflows
    
    Backward Compatibility:
    - DAG workflows unchanged: If no EnhancedSimulation instances exist, behavior matches base class (see run_project in nexus.py)
    - Base class methods called via super() to preserve original functionality
    - New features opt-in: Only EnhancedSimulation instances trigger enhanced functionalities
    """
    
    def __init__(self):
        super().__init__()
        self.filtered_sims: List[Simulation] = []  # Track simulations filtered out
        self._blocked_report: Dict[int, Dict[str, Any]] = {}
        self._blocked_ids: Set[int] = set()
    
    def add_simulations(self, *simulations):
        """Override to filter by definition-time conditionals and blocked dependencies."""
        if len(simulations) == 0:
            simulations = Simulation.all_sims
        if len(simulations) > 0 and not isinstance(simulations[0], Simulation):
            simulations = simulations[0]
        
        # FIRST: Check and mark all simulations as fake if they're machine incompatible
        # This must happen before any other processing
        # Use direct comparison to avoid any issues with is_machine_compatible()
        for sim in simulations:
            if isinstance(sim, EnhancedSimulation):
                required = getattr(sim, 'required_machine', None)
                allowed = getattr(sim, 'allowed_machines', None)
                
                if required or (allowed and len(allowed) > 0):
                    # Get current machine
                    machine_name = None
                    if hasattr(sim, 'job') and sim.job is not None:
                        machine_name = getattr(sim.job, 'machine', None)
                    if machine_name is None:
                        from machines import Job
                        machine_name = Job.machine
                    if machine_name is None:
                        from project_manager import ProjectManager
                        if ProjectManager.machine is not None:
                            machine_name = getattr(ProjectManager.machine, 'name', None)
                    
                    # Direct comparison
                    if required:
                        if machine_name is None or str(machine_name).strip().lower() != str(required).strip().lower():
                            sim.fake_sim = True
                            continue
                    if allowed and len(allowed) > 0:
                        if machine_name is None:
                            sim.fake_sim = True
                            continue
                        machine_normalized = str(machine_name).strip().lower()
                        allowed_normalized = {str(a).strip().lower() for a in allowed}
                        if machine_normalized not in allowed_normalized:
                            sim.fake_sim = True
                            continue
        
        # Filter by definition-time conditionals and machine compatibility
        filtered_sims = []
        self.filtered_sims = []  # Reset filtered list
        blocked_ids = set()  # Track sim_ids of blocked simulations
        
        # First pass: identify directly blocked simulations
        for sim in simulations:
            # Check fake status first - use multiple methods to be sure
            fake_status = sim.fake() if hasattr(sim, 'fake') else False
            fake_attr = getattr(sim, 'fake_sim', False)
            
            # Also check machine compatibility directly (in case fake_sim wasn't set)
            if isinstance(sim, EnhancedSimulation):
                required = getattr(sim, 'required_machine', None)
                if required:
                    machine_name = None
                    if hasattr(sim, 'job') and sim.job is not None:
                        machine_name = getattr(sim.job, 'machine', None)
                    if machine_name is None:
                        from machines import Job
                        machine_name = Job.machine
                    if machine_name is None:
                        from project_manager import ProjectManager
                        if ProjectManager.machine is not None:
                            machine_name = getattr(ProjectManager.machine, 'name', None)
                    if machine_name and str(machine_name).strip().lower() != str(required).strip().lower():
                        fake_attr = True  # Force fake if machine doesn't match
            
            if fake_status or fake_attr:
                # Ensure it's marked as fake
                sim.fake_sim = True
                self.filtered_sims.append(sim)
                blocked_ids.add(sim.simid)
                continue
            
            # Check if it's a EnhancedSimulation with definition-time conditionals
            if isinstance(sim, EnhancedSimulation):
                # Machine filtering - check compatibility
                # This is critical for make_enhanced instances that bypass __init__
                if not sim.is_machine_compatible():
                    sim.fake_sim = True  # Ensure it's marked as fake
                    self.filtered_sims.append(sim)
                    blocked_ids.add(sim.simid)
                    continue
                
                # Definition-time conditional check
                if sim.conditional_enabled:
                    cls = sim.__class__
                    if hasattr(cls, '_conditional_execution_time'):
                        if not cls._conditional_execution_time:  # Definition-time
                            if not sim.evaluate_condition():
                                self.filtered_sims.append(sim)
                                blocked_ids.add(sim.simid)
                                continue
            
            filtered_sims.append(sim)
        # Second pass: filter out simulations that depend on blocked simulations
        # Iterate until no new blocked simulations are found (transitive closure)
        # Create lookup maps for blocked simulations
        blocked_by_id = {}  # simid -> sim
        blocked_by_identifier = {}  # identifier -> sim
        for sim in self.filtered_sims:
            blocked_by_id[sim.simid] = sim
            blocked_by_identifier[sim.identifier] = sim
        
        # Also check all simulations to see if any dependency would be blocked
        # (dependencies might reference simulations not in the input list)
        all_sims_check = Simulation.all_sims if hasattr(Simulation, 'all_sims') else []
        for sim in all_sims_check:
            if sim.simid in blocked_ids:
                blocked_by_id[sim.simid] = sim
                blocked_by_identifier[sim.identifier] = sim
            elif isinstance(sim, EnhancedSimulation):
                # Check if this simulation would be blocked
                if not sim.is_machine_compatible():
                    blocked_by_id[sim.simid] = sim
                    blocked_by_identifier[sim.identifier] = sim
                    if sim not in self.filtered_sims:
                        self.filtered_sims.append(sim)
                        blocked_ids.add(sim.simid)
        
        changed = True
        while changed:
            changed = False
            remaining_sims = []
            for sim in filtered_sims:
                # Check if any dependency is blocked
                has_blocked_dep = False
                blocked_dep_identifier = None
                
                # Check dependency sim objects directly (most reliable)
                for dep in sim.dependencies.values():
                    if hasattr(dep, 'sim'):
                        dep_sim = dep.sim
                        # Check by simid
                        if dep_sim.simid in blocked_ids or dep_sim.simid in blocked_by_id:
                            has_blocked_dep = True
                            blocked_dep_identifier = dep_sim.identifier
                            break
                        # Check by identifier
                        if dep_sim.identifier in blocked_by_identifier:
                            has_blocked_dep = True
                            blocked_dep_identifier = dep_sim.identifier
                            break
                        # Check if dependency sim is in filtered list
                        if dep_sim in self.filtered_sims:
                            has_blocked_dep = True
                            blocked_dep_identifier = dep_sim.identifier
                            break
                
                # Also check dependency_ids set
                if not has_blocked_dep:
                    for dep_id in sim.dependency_ids:
                        if dep_id in blocked_ids or dep_id in blocked_by_id:
                            has_blocked_dep = True
                            # Find the identifier for reporting
                            if dep_id in blocked_by_id:
                                blocked_dep_identifier = blocked_by_id[dep_id].identifier
                            break
                
                # Also check dependencies dict keys
                if not has_blocked_dep:
                    for dep_id in sim.dependencies.keys():
                        if dep_id in blocked_ids or dep_id in blocked_by_id:
                            has_blocked_dep = True
                            if dep_id in blocked_by_id:
                                blocked_dep_identifier = blocked_by_id[dep_id].identifier
                            break
                
                if has_blocked_dep:
                    self.filtered_sims.append(sim)
                    blocked_ids.add(sim.simid)
                    blocked_by_id[sim.simid] = sim
                    blocked_by_identifier[sim.identifier] = sim
                    changed = True
                else:
                    remaining_sims.append(sim)
            
            filtered_sims = remaining_sims
        
        # Final safety check: remove any blocked simulations from filtered_sims
        # (in case they slipped through)
        final_filtered = []
        for sim in filtered_sims:
            if sim.simid not in blocked_ids and not sim.fake() and not getattr(sim, 'fake_sim', False):
                # Double-check machine compatibility for EnhancedSimulations
                if isinstance(sim, EnhancedSimulation):
                    if sim.is_machine_compatible():
                        final_filtered.append(sim)
                    else:
                        # Shouldn't happen, but be safe
                        self.filtered_sims.append(sim)
                        blocked_ids.add(sim.simid)
                else:
                    final_filtered.append(sim)
            else:
                # Shouldn't happen, but be safe
                if sim not in self.filtered_sims:
                    self.filtered_sims.append(sim)
                    blocked_ids.add(sim.simid)
        
        # Final verification: ensure no fake/blocked simulations in final_filtered
        # This is a critical safety check
        verified_filtered = []
        for sim in final_filtered:
            # Triple-check: fake status, blocked_ids, and machine compatibility
            if (sim.simid not in blocked_ids and 
                not sim.fake() and 
                not getattr(sim, 'fake_sim', False)):
                if isinstance(sim, EnhancedSimulation):
                    # One more machine compatibility check
                    if sim.is_machine_compatible():
                        verified_filtered.append(sim)
                    else:
                        # Shouldn't happen, but be safe
                        if sim not in self.filtered_sims:
                            self.filtered_sims.append(sim)
                            blocked_ids.add(sim.simid)
                else:
                    verified_filtered.append(sim)
            else:
                # Shouldn't happen, but be safe
                if sim not in self.filtered_sims:
                    self.filtered_sims.append(sim)
                    blocked_ids.add(sim.simid)
        
        # Call parent method with verified filtered simulations
        # But first, double-check that all simulations are not fake
        final_verified = []
        for sim in verified_filtered:
            # CRITICAL: Check fake status one more time before adding
            if sim.fake() or getattr(sim, 'fake_sim', False):
                # Should never happen, but skip if fake
                continue
            if isinstance(sim, EnhancedSimulation):
                # Double-check machine compatibility
                required = getattr(sim, 'required_machine', None)
                if required:
                    machine_name = None
                    if hasattr(sim, 'job') and sim.job is not None:
                        machine_name = getattr(sim.job, 'machine', None)
                    if machine_name is None:
                        from machines import Job
                        machine_name = Job.machine
                    if machine_name is None:
                        from project_manager import ProjectManager
                        if ProjectManager.machine is not None:
                            machine_name = getattr(ProjectManager.machine, 'name', None)
                    if machine_name and str(machine_name).strip().lower() != str(required).strip().lower():
                        # Machine doesn't match - skip
                        continue
            # Explicitly ensure fake_sim is False
            sim.fake_sim = False
            final_verified.append(sim)
        
        # Snapshot blocked sims (includes downstream dependents) before detaching
        _, snapshot_blocked_ids = self._snapshot_blocked()
        
        # Eliminate blocked simulations from dependency graph before adding new ones
        self._detach_blocked_simulations(snapshot_blocked_ids)
        
        # Final check: ensure parent method won't add any fake simulations
        # The parent checks sim.fake(), so we need to make sure fake_sim is False
        super().add_simulations(final_verified)
        
        # Remove any blocked simulations that might have been added
        # (e.g., if parent method added them despite our filtering)
        self._prune_blocked_from_tracking(snapshot_blocked_ids)
    
    def screen_fake_sims(self):
        """Override to allow blocked simulations (they're intentionally fake)."""
        # For EnhancedProjectManager, blocked simulations are intentionally marked as fake
        # and should not cause an error. We've already filtered them out in add_simulations
        # and run_project, so if any remain, they're legitimate fake sims that should error.
        # But we need to check only simulations that are actually in cascades.
        def collect_fake(sim, fake):
            # Only collect if it's fake AND not intentionally blocked
            if sim.fake():
                # Check if it's a blocked EnhancedSimulation
                if isinstance(sim, EnhancedSimulation):
                    # If it's in our filtered list, it's intentionally blocked - skip
                    if sim in getattr(self, 'filtered_sims', []):
                        return
                fake.append(sim)
        #end def collect_fake
        fake = []
        self.traverse_cascades(collect_fake, fake)
        if len(fake) > 0:
            msg = 'fake/temporary simulation objects detected in cascade\nthis is a developer error\nlist of fake sims and directories:\n'
            for sim in fake:
                msg += '  {0:>8}  {1}\n'.format(sim.simid, sim.locdir)
            #end for
            self.error(msg)
        #end if
    #end def screen_fake_sims
    
    def run_project(self, status=False, status_only=False):
        """Override to ensure blocked simulations are removed before running."""
        # Before running, ensure all blocked simulations are removed
        # This is critical for generate_only and other operations
        _, blocked_ids = self._snapshot_blocked()
        self._detach_blocked_simulations(blocked_ids)
        self._prune_blocked_from_tracking(blocked_ids)
        # Now call parent method
        super().run_project(status=status, status_only=status_only)
    
    def init_cascades(self):
        """Override to ensure blocked simulations are removed before initializing cascades."""
        # Remove blocked simulations BEFORE initializing cascades
        # This is critical because cascades are created in add_simulations
        _, blocked_ids = self._snapshot_blocked()
        self._detach_blocked_simulations(blocked_ids)
        self._prune_blocked_from_tracking(blocked_ids)
        # Now call parent method
        super().init_cascades()
    
    def traverse_cascades(self, operation=None, *args, **kwargs):
        """Override to use standard DAG traversal."""
        # Reset wait_ids for all cascades
        for cascade in self.cascades.values():
            cascade.reset_wait_ids()
        
        # Use standard DAG traversal (all workflows are DAGs)
        for cascade in self.cascades.values():
            cascade.traverse_cascade(operation or trivial, *args, **kwargs)
        
    def _detach_blocked_simulations(self, blocked_ids: Set[int]):
        """Remove blocked simulations from dependency graph (dependents and dependencies)."""
        if not blocked_ids:
            return
        
        sim_lookup = {}
        all_sims = getattr(Simulation, 'all_sims', [])
        for sim in all_sims:
            if hasattr(sim, 'simid'):
                sim_lookup[sim.simid] = sim
        
        for sim_id in blocked_ids:
            blocked_sim = sim_lookup.get(sim_id)
            if blocked_sim is None:
                continue
            
            if hasattr(blocked_sim, 'eliminate'):
                blocked_sim.eliminate()
            else:
                # Manual removal as fallback
                dependencies = list(getattr(blocked_sim, 'dependencies', {}).values())
                for dep in dependencies:
                    upstream = getattr(dep, 'sim', None)
                    if upstream and hasattr(upstream, 'dependents'):
                        upstream.dependents.pop(sim_id, None)
                dependents = list(getattr(blocked_sim, 'dependents', {}).values())
                for dependent in dependents:
                    if hasattr(dependent, 'dependencies') and sim_id in dependent.dependencies:
                        dependent.undo_depends(blocked_sim)
            
            # Ensure blocked simulations remain blocked/fake
            blocked_sim.block = True
            blocked_sim.block_subcascade = True
            blocked_sim.fake_sim = True
    
    def _prune_blocked_from_tracking(self, blocked_ids: Set[int]):
        """Remove blocked simulations from manager tracking dictionaries."""
        if not blocked_ids:
            return
        
        for sim_id in list(self.simulations.keys()):
            sim = self.simulations[sim_id]
            if (sim_id in blocked_ids or
                sim.fake() or
                getattr(sim, 'fake_sim', False) or
                (isinstance(sim, EnhancedSimulation) and not sim.is_machine_compatible())):
                del self.simulations[sim_id]
        
        for sim_id in list(self.cascades.keys()):
            if sim_id in blocked_ids:
                del self.cascades[sim_id]
        
        for sim_id in list(self.progressing_cascades.keys()):
            if sim_id in blocked_ids:
                del self.progressing_cascades[sim_id]
    
    def _snapshot_blocked(self):
        """Capture blocked simulations before mutating the workflow graph."""
        blocked_sims, blocked_ids = self._get_blocked_simulations()
        for sim, reasons in blocked_sims:
            entry = dict(
                simid=sim.simid,
                identifier=str(sim.identifier),
                path=str(sim.path),
                reasons=list(reasons)
            )
            self._blocked_report[sim.simid] = entry
        self._blocked_ids.update(blocked_ids)
        return blocked_sims, blocked_ids
    
    def _get_blocked_simulations(self):
        """Identify all blocked simulations (directly and by dependencies)."""
        blocked_sims = []
        directly_blocked = set()  # Track sim_ids of directly blocked simulations
        
        # Check simulations that were filtered out during add_simulations
        sims_to_check = []
        sims_to_check.extend(getattr(self, 'filtered_sims', []))
        
        # Also check all_sims for any EnhancedSimulations that might be blocked
        all_sims = Simulation.all_sims if hasattr(Simulation, 'all_sims') else []
        # Create a set of simids to avoid recursion issues with 'in' operator
        sims_to_check_simids = {s.simid for s in sims_to_check if hasattr(s, 'simid')}
        for sim in all_sims:
            if isinstance(sim, EnhancedSimulation):
                # Use simid comparison to avoid recursion issues with __eq__
                if hasattr(sim, 'simid') and sim.simid not in sims_to_check_simids:
                    # Check if it would be blocked
                    if not sim.is_machine_compatible() or (sim.conditional_enabled and 
                        hasattr(sim.__class__, '_conditional_execution_time') and
                        not sim.__class__._conditional_execution_time and
                        not sim.evaluate_condition()):
                        sims_to_check.append(sim)
                        sims_to_check_simids.add(sim.simid)
        
        # First pass: identify directly blocked simulations
        for sim in sims_to_check:
            reasons = []
            
            if isinstance(sim, EnhancedSimulation):
                # Check machine compatibility
                if not sim.is_machine_compatible():
                    machine_name = getattr(sim.job, 'machine', None) if hasattr(sim, 'job') and sim.job else 'None'
                    required = getattr(sim, 'required_machine', None)
                    allowed = getattr(sim, 'allowed_machines', None)
                    
                    if required:
                        reasons.append(f"required_machine='{required}' (current: '{machine_name}')")
                    elif allowed and len(allowed) > 0:
                        reasons.append(f"allowed_machines={allowed} (current: '{machine_name}')")
                    else:
                        reasons.append(f"machine incompatible (current: '{machine_name}')")
                
                # Check definition-time conditionals
                if sim.conditional_enabled:
                    cls = sim.__class__
                    if hasattr(cls, '_conditional_execution_time'):
                        if not cls._conditional_execution_time:  # Definition-time
                            try:
                                if not sim.evaluate_condition():
                                    reasons.append("definition-time conditional failed")
                            except Exception as e:
                                reasons.append(f"definition-time conditional error: {str(e)}")
                
                # Check if skipped
                if getattr(sim, 'skipped', False):
                    reasons.append("skipped")
            
            if reasons:
                blocked_sims.append((sim, reasons))
                directly_blocked.add(sim.simid)
        
        # Second pass: check for simulations blocked by dependencies
        # Iterate until no new blocked simulations are found (transitive closure)
        changed = True
        while changed:
            changed = False
            for sim in all_sims:
                if sim.simid in directly_blocked:
                    continue
                
                # Check if any dependency is blocked
                blocked_deps = []
                for dep_id, dep in sim.dependencies.items():
                    if dep_id in directly_blocked:
                        dep_sim = dep.sim
                        blocked_deps.append(dep_sim.identifier)
                
                if blocked_deps:
                    reasons = [f"dependency blocked: {', '.join(blocked_deps)}"]
                    blocked_sims.append((sim, reasons))
                    directly_blocked.add(sim.simid)
                    changed = True
        
        return blocked_sims, directly_blocked
    
    def write_simulation_status(self):
        """Override to filter blocked simulations from normal status and add blocked section."""
        # Snapshot current state to ensure any new blocked sims are captured
        self._snapshot_blocked()
        blocked_ids = set(self._blocked_ids)
        
        # Temporarily remove blocked simulations from self.simulations for status output
        original_simulations = self.simulations.copy()
        filtered_simulations = obj()
        for sim_id, sim in self.simulations.items():
            if sim_id not in blocked_ids:
                filtered_simulations[sim_id] = sim
        
        # Temporarily replace simulations dict
        self.simulations = filtered_simulations
        
        # Call parent method (will only show non-blocked simulations)
        super().write_simulation_status()
        
        # Restore original simulations dict
        self.simulations = original_simulations
        
        # Report blocked simulations
        if self._blocked_report:
            blocked_entries = sorted(self._blocked_report.values(), key=lambda x: x['identifier'])
            self.log('\n==== blocked simulations (will not run) ====', n=1)
            for entry in blocked_entries:
                reason_str = '; '.join(entry['reasons'])
                self.log(f"  {entry['identifier']:30s} ({entry['path']:40s}) - {reason_str}", n=2)
            self.log('', n=1)
    
    def status_line(self, sim, extra=''):
        """Override to show iteration_count for simulations with loops."""
        from enhanced_simulation import EnhancedSimulation
        indicators = ('setup', 'sent_files', 'submitted', 'finished', 'got_output', 'analyzed')
        stats = sim.tuple(*indicators)
        status = ''
        for stat in stats:
            status += str(int(stat))
        if sim.process_id is None:
            pid = '------'
        else:
            pid = sim.process_id
        # Add enhanced information for EnhancedSimulation instances
        if isinstance(sim, EnhancedSimulation):
            info_parts = []
            
            # Loop information
            if sim.loop_enabled:
                iter_info = f"iter={sim.iteration_count}"
                if sim.loop_max_iterations is not None:
                    iter_info += f"/{sim.loop_max_iterations}"
                info_parts.append(iter_info)
            
            # Conditional information
            cond_info = sim.get_conditional_info()
            if cond_info:
                info_parts.append(cond_info)
            
            # Merge strategy information
            if sim.merge_selector is not None:
                merge_info = f"merge={sim.merge_strategy}"
                selector_name = getattr(sim.merge_selector, '__name__', 'custom')
                merge_info += f":{selector_name}"
                info_parts.append(merge_info)
            
            # Build status line with enhanced info
            if info_parts:
                enhanced_info = ' '.join(info_parts)
                sline = '{0}  {1}  {2:<8}  {3:<6}  {4} [{5}]'.format(
                    status, str(int(sim.failed)), pid, sim.identifier, sim.locdir, enhanced_info
                )
            else:
                sline = '{0}  {1}  {2:<8}  {3:<6}  {4}'.format(
                    status, str(int(sim.failed)), pid, sim.identifier, sim.locdir
                )
        else:
            sline = '{0}  {1}  {2:<8}  {3:<6}  {4}'.format(
                status, str(int(sim.failed)), pid, sim.identifier, sim.locdir
            )
        self.log(sline, extra, n=2)

