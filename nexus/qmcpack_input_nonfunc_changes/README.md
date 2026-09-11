# Deferred Nexus `qmcpack_input.py` changes

This directory intentionally keeps the Nexus-side work removed from PR #6184.
The active `nexus/qmcpack_input.py` and
`nexus/tests/test_qmcpack_input.py` have been restored to `develop`, so the
PR is limited to QMCPACK's C++ removal of non-functional inputs.

## Snapshots

- `qmcpack_input.py` is the PR-version snapshot, including the resolution of
  the `schema_e` alias-handling merge conflict.
- `qmcpack_input.py.develop` is the exact `develop` version restored to the
  active Nexus tree (base commit `ef1dd3bbc627f82914da75bd21d1775c4ec50e9c`).

These are deliberately adjacent so the deferred change can be inspected with:

```bash
diff -u qmcpack_input_nonfunc_changes/qmcpack_input.py.develop \
        qmcpack_input_nonfunc_changes/qmcpack_input.py
```

At archival time, that direct comparison is 1,346 additions and 400 deletions.

## What the deferred change does

1. Expands the declarative Nexus XML schema: accepted elements, attributes,
   parameters, aliases, plural collections, and newly represented XML nodes.
2. Makes `QIxml.init_from_xml` resolve element aliases through `schema_e`
   before looking up types, plural collections, or permitted elements. This
   is the behavior retained when resolving the PR's merge conflict.
3. Adds a broad Nexus input-contract test,
   `test_qixml_live_and_unsupported_parameters`, in the removed companion
   change to `nexus/tests/test_qmcpack_input.py`. It distinguishes recognized
   live inputs from retired/unsupported ones.

## Roll-forward / rollback notes

- Do not copy this snapshot over the active file blindly. Rebase the two-file
  comparison above onto the then-current `develop` version and resolve schema
  additions against current parser behavior.
- The schema lists duplicate knowledge held by QMCPACK's C++ XML parsers.
  Treat this as a compatibility-contract update and validate it with the
  Nexus input tests when reintroducing it.
- If legacy input is blocked before this work returns, first reassess the
  retired-input cases in the removed test: tests that only prove legacy
  no-ops may no longer justify permanent maintenance.
- The original active PR snapshot before this archival split is merge commit
  `93f0f289400848330e9b6f67dc5947f0227f7931`; it is useful as an additional
  historical reference if the local snapshots need cross-checking.
