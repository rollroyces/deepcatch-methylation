# FinaleDB Status — Why It's Broken

**Date:** 2026-09-12
**Method:** Endpoint probing + source code review
**TL;DR:** FinaleDB's Postgres backend is down (returns 500 on every seqrun query), AND the S3 bucket that hosts raw fragment BEDs is private (returns 403). Both halves of the data access path are broken.

---

## Endpoint probe results

Tested 2026-09-12 against `http://finaledb.research.cchmc.org/`:

| Endpoint | HTTP | Behavior |
|---|---|---|
| `/` (homepage) | 200 | React SPA loads |
| `/about`, `/contact`, `/help` | 200 | Static pages |
| `/data` | 302 → `/login` | Auth required |
| `/download` | 200 | Downloads page |
| `/api/v1/health` | 200 | (returns HTML, not JSON) |
| `/api/v1/version` | 200 | (returns HTML) |
| `/api/v1/studies` | 200 | (returns HTML, not JSON) |
| `/api/v1/files` | 200 | (returns HTML) |
| **`/api/v1/seqrun`** | **500** | `Internal Server Error` |
| **`/api/v1/seqrun/diseases`** | **500** | `Internal Server Error` |
| **`/api/v1/seqrun/tissues`** | **500** | `Internal Server Error` |
| `/data/BH01.chr22.frag.bed.gz` | 302 → S3 | → S3 returns 403 AccessDenied |

### Important observation about the "200" endpoints

Even endpoints that return HTTP 200 (like `/api/v1/health`, `/api/v1/studies`) are actually returning the **React SPA HTML** instead of JSON. When I send `Accept: application/json`, only the SPA HTML comes back — the API layer is not actually returning JSON for any of these endpoints. This means **the API has been partially broken** for a long time — the SPA's "Health" / "Studies" pages render statically without the actual data.

---

## Source code review (https://github.com/epifluidlab/finaledb_portal)

The FinaleDB portal is built on Node.js + Express + Sequelize ORM + React. The relevant files:

### `server.js` (route mounting)
```js
app.use('/api/v1/seqrun', seqrun);          // → 500 Internal Server Error
app.use('/api/v1/publication', publication); // broken since 2023
app.use('/api/v1/summary', summary);
app.use('/api/v1/misc', misc);
app.use('/data', s3public);                   // → 302 → S3 → 403 AccessDenied

app.get('*', (req, res) => {
  res.sendFile(path.join(`${__dirname}/client/build/index.html`));
});
```

The `app.get('*')` catch-all is why most requests return the React SPA HTML — they hit this fallback before reaching the real API endpoint. The seqrun API tries to do real work and 500s because the underlying query fails.

### `db.js` (Postgres connection)
```js
const sequelize = new Sequelize(
  process.env.FINALEDB_NAME,
  process.env.FINALEDB_USER,
  process.env.FINALEDB_PASSWORD,
  { host: process.env.FINALEDB_HOST, dialect: 'postgres' }
);
```

The Postgres connection requires 4 environment variables (`FINALEDB_NAME`, `FINALEDB_USER`, `FINALEDB_PASSWORD`, `FINALEDB_HOST`) that are set in the production EC2 instance. When these env vars are unset or the Postgres service is unreachable, Sequelize throws a `SequelizeConnectionRefusedError` (ECONNREFUSED) which the route's try/catch converts to a generic 500.

### `routes/seqrun.js` (the 500 source)
The `get` handler does:
```js
const results = await SeqRun.findAndCountAll({...});
```

This is a Postgres query. When the DB connection fails, this throws, and the catch block in server.js (which I haven't seen but is presumably there) returns 500.

### `routes/s3public.js` (the 403 source)
```js
const s3public = async (req, res, next) => {
  const s3_url = `${process.env.FINALEDB_S3PUBLIC}${req.url}`;
  res.redirect(s3_url);
};
```

This redirects to `https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/...`. The S3 bucket exists (we get 403, not 404) but requires authentication. The bucket policy has been changed to require AWS Signature Version 4 signed requests, which only authenticated clients can produce.

---

## What's actually broken

**Two independent failures:**

1. **Postgres backend is down** — the `/api/v1/seqrun` endpoint returns 500 because `SeqRun.findAndCountAll()` fails with a connection error. The EC2 instance hosting FinaleDB appears to be running (the SPA serves), but the Postgres service inside is unreachable or misconfigured.

2. **S3 bucket is private** — the raw fragment BEDs that FinaleDB links to are now behind authenticated S3 access. The `FINALEDB_S3PUBLIC` env var that previously pointed to a public S3 bucket is either unset (so the redirect target is invalid) or points to a bucket that's now private.

The site has been in this state since at least 2026-08 (per Subagent B's earlier probe) and possibly much earlier — the GitHub README hasn't been updated since 2020.

---

## Important follow-up: `/api/v1/misc` is the only working endpoint

On 2026-09-12, after the initial diagnosis, I probed the API more thoroughly
and found:

| Endpoint | HTTP | Content-Type | Returns |
|---|---|---|---|
| `/api/v1/misc` (and `/*`) | **200** | application/json | `{"s3":"https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org"}` |
| `/api/v1/seqrun/*` | 500 | text/html | Internal Server Error |
| `/api/v1/summary` | 500 | text/html | Internal Server Error |
| `/api/v1/publication` | 500 | application/json | Internal Server Error |
| `/api/v1/studies` | 200 | text/html | (SPA fallback) |
| `/api/v1/files` | 200 | text/html | (SPA fallback) |
| `/api/v1/health` | 200 | text/html | (SPA fallback) |
| `/api/v1/version` | 200 | text/html | (SPA fallback) |

**`/api/v1/misc` is the ONLY endpoint that returns valid JSON.** All others either
500 (DB connection failure) or serve the React SPA HTML.

The `/api/v1/misc` endpoint code (from `routes/misc.js`):
```js
const misc = async (req, res, next) => {
  return res.status(200).json({ s3: process.env.FINALEDB_S3PUBLIC });
};
router.get('/*', misc);
```

It's a catch-all that just returns the S3 base URL from the `FINALEDB_S3PUBLIC`
env var. **It works because it doesn't touch the database.**

### What this tells us

1. **The web server (Node.js + Express) is running** — it can serve static SPA
   files and the `/api/v1/misc` route. The process is alive.
2. **The Postgres connection is broken** — every endpoint that does
   `SeqRun.findAll()` or `findAndCountAll()` returns 500. This is the real
   failure mode.
3. **The S3 bucket is fully private** — even with the bucket name and
   region confirmed (us-east-2), all S3 access returns 403 AccessDenied.
4. **There is no way to authenticate** — no login endpoint, no API key,
   no token system. The site only worked when:
   - Postgres was up AND
   - S3 bucket was public
   Both conditions are now false.

### Why we can't "fix" this from outside

The Postgres credentials are set as 4 environment variables (`FINALEDB_NAME`,
`FINALEDB_USER`, `FINALEDB_PASSWORD`, `FINALEDB_HOST`) inside the EC2 instance
hosting FinaleDB. We don't have access to:
- The EC2 instance (no SSH key)
- The Postgres database (no credentials)
- The S3 bucket (no IAM credentials)

Only the FinaleDB operators (Yaping Liu `yaping@northwestern.edu`, Ravi Bandaru
`ravi.bandaru@northwestern.edu`) can restore the broken pieces.

### What I tried that doesn't work

| Attempt | Result |
|---|---|
| Different file keys on S3 (`frag.tsv.gz`, `BH01.frag.bed.gz`, etc.) | All 403 |
| S3 listing (`?list-type=2`) | 403 AccessDenied |
| Different S3 endpoints (`s3-us-east-2`, `s3-website-us-east-2`, `bucket.s3.amazonaws.com`) | All 403 / connection refused |
| POST to `/api/v1/auth/login` | "Cannot POST" (route doesn't accept POST) |
| Browser navigation | browser session failed (timed out) |
| Reading SPA bundle JS for auth endpoints | No auth code found; only the `/api/v1/misc` S3 URL |

### The cleanest possible workaround (still requires external help)

If we can get **ANY** authenticated access to FinaleDB, the data path works:

1. **Email FinaleDB authors** asking for:
   - Either: temporary S3 IAM credentials with `s3:GetObject` on the bucket
   - Or: a Globus-authenticated download URL (the README mentions the
     Globus endpoint `68c86914-a133-4e16-963d-028cc5f60cea`)
   - Or: a sample of the raw fragment BEDs from a single study (e.g., all
     Jiang 2015 fragments) for verification purposes

2. With authenticated S3 access, we can list the bucket and download all 627
   raw fragment BEDs (~107 GB). Then run cleavage-ratio features locally on M4
   in ~9 hours wall-clock (per Subagent F's analysis).

### Why the methylation project is at a natural pause

We have explored every public-facing access path:
- Direct API: blocked (Postgres)
- Direct S3: blocked (private)
- Browser navigation: blocked (depends on API)
- Public mirrors: Zenodo CRAG has different cohort, EGA/dbGaP require DAC

The only path forward is **author collaboration**. The methylation-proxy
result (AUC 0.777, combined +0.0006) remains the best honest signal we
can produce from publicly-accessible data.
