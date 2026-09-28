"""Shared zip-extraction safety. Single home of the member-path traversal guard
(previously a documented local copy in raw_antismash_triage; see that module's
docstring for the mamey-must-not-import-tools rationale, which this placement keeps)."""
from __future__ import annotations
import os
import stat
import unicodedata
import zipfile
from pathlib import Path


def is_regular_file_member(member: zipfile.ZipInfo) -> bool:
    """True only for a ZIP member that can safely be consumed as file data.

    ZIPs made on DOS/Windows and some Python fixtures do not carry Unix file-type
    bits; a zero type is therefore treated as an ordinary file.  Explicit Unix
    special types (symlink, device, FIFO, socket) are not file data and must not
    be admitted merely because their names end in ``.gbk`` or ``.json``.
    """
    if member.is_dir():
        return False
    mode = (member.external_attr >> 16) & 0xFFFF
    file_type = stat.S_IFMT(mode)
    return file_type in (0, stat.S_IFREG)


def is_macos_metadata_name(name: str) -> bool:
    """True for Finder metadata entries that must not be consumed as data."""
    base = name.rsplit("/", 1)[-1]
    return (
        name.startswith("__MACOSX/")
        or "/__MACOSX/" in name
        or base.startswith("._")
        or base == ".DS_Store"
    )


def regular_file_names(z: zipfile.ZipFile) -> list[str]:
    """Return the effective regular data-file namespace of *z*.

    ``ZipFile.open(name)`` and ``getinfo(name)`` resolve duplicate names to the
    last central-directory entry.  Mirror that rule before filtering so an
    earlier regular entry cannot make a later special entry with the same name
    look admissible, and so parsers do not count duplicate names twice. Finder
    resource-fork trees, AppleDouble sidecars and ``.DS_Store`` are regular
    files at the ZIP layer but metadata rather than antiSMASH input data.
    """
    effective: dict[str, zipfile.ZipInfo] = {}
    for member in z.infolist():
        effective[member.filename] = member
    return [name for name, member in effective.items()
            if is_regular_file_member(member)
            and not is_macos_metadata_name(name)]


def nonregular_file_names(z: zipfile.ZipFile) -> list[str]:
    """Names of explicit non-directory special members, in archive order."""
    return [member.filename for member in z.infolist()
            if not member.is_dir() and not is_regular_file_member(member)]


def duplicate_member_names(z: zipfile.ZipFile) -> list[str]:
    """Sorted member names that occur more than once in the central directory."""
    seen: set[str] = set()
    duplicates: set[str] = set()
    for member in z.infolist():
        if member.filename in seen:
            duplicates.add(member.filename)
        seen.add(member.filename)
    return sorted(duplicates)


def _zipfile_extraction_target(dest_abs: Path, member: zipfile.ZipInfo) -> Path:
    """Return the target path produced by ``ZipFile.extract`` normalization.

    The standard-library extractor drops empty, ``.`` and ``..`` components
    instead of resolving ``..`` as a parent traversal.  Collision preflight
    must mirror that behavior or two distinct member strings can overwrite the
    same target even though their ordinary ``Path.resolve`` results differ.
    """
    arcname = member.filename
    if os.path.sep != "/":
        arcname = arcname.replace("/", os.path.sep)
    if os.path.altsep and os.path.altsep != "/":
        arcname = arcname.replace(os.path.altsep, os.path.sep)
    _, arcname = os.path.splitdrive(arcname)
    invalid_parts = ("", os.path.curdir, os.path.pardir)
    arcname = os.path.sep.join(
        part for part in arcname.split(os.path.sep) if part not in invalid_parts
    )
    if os.path.sep == "\\":
        arcname = zipfile.ZipFile._sanitize_windows_name(arcname, os.path.sep)
    if not arcname and not member.is_dir():
        raise ValueError(f"unsafe empty zip extraction target: {member.filename}")
    return Path(os.path.normpath(os.path.join(dest_abs, arcname))).resolve()


def _portable_extraction_key(dest_abs: Path, target: Path) -> tuple[str, ...]:
    """Conservative cross-filesystem key for an admitted extraction target.

    NFKC plus Unicode casefold deliberately rejects some distinct names that a
    case-sensitive filesystem could store.  That stable over-rejection avoids
    silent overwrites when the same archive is moved to common case-insensitive
    or Unicode-normalizing filesystems.  It is not an emulator for every
    filesystem-specific alias rule.
    """
    relative = target.relative_to(dest_abs)
    return tuple(
        unicodedata.normalize("NFKC", part).casefold()
        for part in relative.parts
    )


def safe_extract_all(z: zipfile.ZipFile, dest: Path | str) -> None:
    """Extract ordinary files/directories after an all-members safety preflight.

    Validation completes before the first write.  Absolute/traversal paths,
    duplicate names, and Unix special members are refused rather than partially
    extracting an archive whose remaining members are unsafe or ambiguous.
    """
    dest_abs = Path(dest).resolve()
    members = z.infolist()
    duplicates = duplicate_member_names(z)
    if duplicates:
        raise ValueError(f"unsafe duplicate zip member path: {duplicates[0]}")
    resolved_members: list[tuple[zipfile.ZipInfo, Path]] = []
    target_owner: dict[Path, zipfile.ZipInfo] = {}
    portable_owner: dict[tuple[str, ...], zipfile.ZipInfo] = {}
    for member in members:
        raw_target = (dest_abs / member.filename).resolve()
        if not (
            raw_target == dest_abs
            or str(raw_target).startswith(str(dest_abs) + os.sep)
        ):
            raise ValueError(f"unsafe zip member path: {member.filename}")
        if not member.is_dir() and not is_regular_file_member(member):
            raise ValueError(f"unsafe non-regular zip member: {member.filename}")
        target = _zipfile_extraction_target(dest_abs, member)
        if not (target == dest_abs or str(target).startswith(str(dest_abs) + os.sep)):
            raise ValueError(f"unsafe zip extraction target: {member.filename}")
        prior = target_owner.get(target)
        if prior is not None:
            raise ValueError(
                "unsafe colliding zip extraction target: "
                f"{prior.filename} and {member.filename}"
            )
        portable_key = _portable_extraction_key(dest_abs, target)
        portable_prior = portable_owner.get(portable_key)
        if portable_prior is not None:
            raise ValueError(
                "unsafe portable zip extraction target collision: "
                f"{portable_prior.filename} and {member.filename}"
            )
        target_owner[target] = member
        portable_owner[portable_key] = member
        resolved_members.append((member, target))

    # A regular-file member cannot also be an ancestor directory needed by a
    # second member.  ZipFile discovers that conflict only while extracting,
    # after it may already have written earlier members.  Reject the complete
    # target shape during preflight so an invalid archive leaves no partial
    # extraction behind.
    file_targets = {target for member, target in resolved_members if not member.is_dir()}
    portable_file_targets = {
        _portable_extraction_key(dest_abs, target)
        for member, target in resolved_members
        if not member.is_dir()
    }
    for member, target in resolved_members:
        parent = target.parent
        while parent != dest_abs and str(parent).startswith(str(dest_abs) + os.sep):
            blocker = target_owner.get(parent)
            if parent in file_targets and blocker is not None:
                raise ValueError(
                    "unsafe zip member file/directory collision: "
                    f"{blocker.filename} blocks {member.filename}"
                )
            parent = parent.parent
        portable_key = _portable_extraction_key(dest_abs, target)
        for width in range(1, len(portable_key)):
            parent_key = portable_key[:width]
            blocker = portable_owner.get(parent_key)
            if parent_key in portable_file_targets and blocker is not None:
                raise ValueError(
                    "unsafe portable zip member file/directory collision: "
                    f"{blocker.filename} blocks {member.filename}"
                )
    for member in members:
        z.extract(member, dest_abs)
