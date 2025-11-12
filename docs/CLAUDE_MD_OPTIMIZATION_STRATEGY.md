# CLAUDE.md Optimization Strategy

## Current State Analysis
- **Size**: 54 lines, 4KB (already compact)
- **Content**: Mix of critical alerts, commands, architecture summary, and recent updates
- **Issues**: Recent updates section will grow over time, potentially affecting performance

## Optimization Goals
1. Keep CLAUDE.md under 40 lines for optimal Claude Code performance
2. Retain only immediately actionable information
3. Move historical and detailed content to appropriate documents
4. Create clear navigation structure

## Content Reorganization Plan

### 1. CLAUDE.md (Keep - Ultra Concise)
**Target**: < 40 lines
- Project identification (1 line)
- Critical alerts (PERV, PMDA, Security) - 4 lines
- Essential commands (test, format, run) - 5 lines
- Quick architecture reference - 2 lines
- Key files reference - 3 lines
- Documentation index table - 5 lines
- Total: ~25-30 lines

### 2. Move to README.md
- Detailed project overview
- Full architecture description
- Setup instructions
- Public-facing documentation

### 3. Move to docs/CHANGELOG.md
- All "Recent Updates" content
- Version history
- Protocol updates
- Session summaries

### 4. Move to docs/claude-sessions/
- Detailed session logs (already exists)
- Protocol change implementations
- Keep only latest session reference in CLAUDE.md if needed

### 5. Create docs/QUICK_REFERENCE.md
- Extended command list
- Common workflows
- Debugging tips
- AWS operations

### 6. Existing Structure (Keep)
- docs/COMMANDS.md - Detailed command reference
- docs/PATTERNS.md - Code patterns and conventions
- docs/ARCHITECTURE.md - Detailed architecture
- docs/PROTOCOLS_GUIDE.md - Sample prep protocols
- docs/DEVELOPMENT_GUIDE.md - Development guide

## Implementation Steps

1. **Extract Recent Updates**
   - Move entire "Recent Updates" section to CHANGELOG.md
   - Add entry with timestamp

2. **Simplify Architecture Section**
   - Keep one-line pipeline summary
   - Move detailed AWS pattern to ARCHITECTURE.md

3. **Streamline Critical Context**
   - Combine warnings into bullet points
   - Remove explanatory text

4. **Create Quick Reference**
   - Extract extended commands
   - Add common troubleshooting

5. **Update Documentation Table**
   - Add link to CHANGELOG.md
   - Add link to QUICK_REFERENCE.md
   - Ensure all links are working

## Expected Results

### Before
- 54 lines, mixed content types
- Growing "Recent Updates" section
- Potential performance impact over time

### After
- ~30 lines, focused on essentials
- Historical data in appropriate files
- Optimized for Claude Code parsing
- Clear navigation to detailed docs

## Validation Checklist
- [ ] CLAUDE.md < 40 lines
- [ ] All critical information retained
- [ ] No broken links
- [ ] Historical data preserved in CHANGELOG
- [ ] README.md updated with project details
- [ ] All documentation accessible via index