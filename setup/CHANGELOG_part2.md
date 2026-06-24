---
title: Part 2 changes
date: 2026-06-24
author: Fred Jaya
---

# Part 2 changes

> Informal notes on what changed and why, for anyone picking this up or delivering it for the first time.

## General

Coming to this as a first-time trainer, changes to the content were made to make the content clearer for delivery, and providing explicit instructions for the exercises.

- **Explicit exercise tasks:** code is provided, or tasks broken down to minimise confusion for learners (and me!)
- **Checkpoint files (`.qs2`):** added after each section so learners can restore a known-good state if they go off-script during exercises. Paired with a catchup script so nobody gets left behind. Can encourage learners to "mess things up" and have a play with a safety net
- **Learning objectives per lesson:** each section now opens with a short list of what learners should be able to do by the end, helped identify the key concepts, trim down content to find more time for self-guided exercises.
- **"See more" folded sections:** deeper technical detail is collapsed by default. The content is still there for curious learners, but doesn't slow down the main flow for the target audience.

## Introduction (`00_intro.Rmd`)

**New page.** Recaps Part 1 content so learners don't have to hold everything in their head across sessions:

- What the `seurat_obj` is and where it came from
- The experimental design and biological hypotheses we'll be testing in Part 2

## Clustering (`02_clustering.Rmd`)

- Reduced the number of clustering resolutions explored — fewer options to compare means learners spend more time understanding the concept and less time waiting for plots to render.

## Cluster Markers (`03_clustermarkers.Rmd`)

- Trimmed visualisation types down to those that directly serve the learning objectives — removed extras that added noise without adding understanding.
- Exercises now walk through the *thinking*, not just the mechanics:
  - Prompt learners to consider which cell types are relevant to their sample
  - Ask them to think about how they'd find canonical markers and why those markers are meaningful
  - Step through manual annotation for a handful of clusters as a worked example

## SingleR (`04_singler.Rmd`)

- Added an explicit comparison of **reference-based annotation** (SingleR) vs. **marker-based annotation** (manual), covering the trade-offs of each approach so learners can make an informed choice for their own data.

## Differential Expression (`05_de2.Rmd`)

- Integrated a **pre-filtering step** before running DEGs to reduce noise and improve results.
- Adopted the **limma-edgeR workflow** as the primary DE method.
- Structured the section around two complementary hypotheses worked through interactively:
  - Global DE (across all cells)
  - Cell-type-specific DE

## Wrapping Up & Next Steps (`06_next_steps.Rmd`)

**New section.** Covers analyses that go beyond the workshop scope but are common next steps in real projects — gives learners a map of where to go after the course ends.
