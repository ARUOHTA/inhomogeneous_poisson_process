---
name: ejdb-open-paper
description: >
  Use this skill whenever the user asks to open, view, or access
  an academic paper, journal article, or DOI.
  This skill enforces the use of Kyoto University EJDB authentication
  via the official working bookmarklet method.
---

# EJDB Academic Paper Access Skill (Kyoto University)

## Purpose

This skill ensures that academic papers are accessed through
Kyoto University’s Electronic Journal & Database (EJDB) system
using the officially working JavaScript bookmarklet method.

The final step (clicking the bookmarklet and adding to Paperpile)
is always performed manually by the user.

## Trigger Conditions

Always use this skill when the user requests actions such as:
- "Open this paper"
- "View this article"
- "Access this DOI"
- "Read the journal paper"
- Any instruction that implies opening an academic paper URL

## Critical Rule

**Do NOT generate custom EJDB URLs.**  
**Do NOT use URL encoding.**

Kyoto University’s current EZproxy system only works reliably
with the following JavaScript bookmarklet pattern.

## Instructions

When this skill is triggered, follow these steps strictly:

1. Identify the target academic URL (DOI or publisher page).
2. Present the original URL clearly.
3. Instruct the user to open the URL normally in the browser.
4. Explicitly instruct the user to click their EJDB bookmarklet
   that uses the following JavaScript code:

javascript:void(location.href='https://kyoto-u.idm.oclc.org/login?url=%27+location.href);


5. Explain that after authentication, the user can add the paper
   to Paperpile with one click.

## Output Format

Use the following format exactly:

### Paper URL
- ORIGINAL_URL

### How to open via EJDB
1. Open the paper URL in your browser.
2. Click your EJDB bookmarklet (Kyoto University).
3. After authentication, add the paper to Paperpile.

## Notes

- Do not attempt to open the browser automatically.
- Do not attempt to download PDFs.
- Do not modify or rewrite the JavaScript code.
- The user performs the final click and Paperpile import manually.

