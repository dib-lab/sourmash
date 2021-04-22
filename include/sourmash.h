/* c bindings to the sourmash library */

#ifndef SOURMASH_H_INCLUDED
#define SOURMASH_H_INCLUDED

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

enum HashFunctions {
  HASH_FUNCTIONS_MURMUR64_DNA = 1,
  HASH_FUNCTIONS_MURMUR64_PROTEIN = 2,
  HASH_FUNCTIONS_MURMUR64_DAYHOFF = 3,
  HASH_FUNCTIONS_MURMUR64_HP = 4,
};
typedef uint32_t HashFunctions;

enum SourmashErrorCode {
  SOURMASH_ERROR_CODE_NO_ERROR = 0,
  SOURMASH_ERROR_CODE_PANIC = 1,
  SOURMASH_ERROR_CODE_INTERNAL = 2,
  SOURMASH_ERROR_CODE_MSG = 3,
  SOURMASH_ERROR_CODE_UNKNOWN = 4,
  SOURMASH_ERROR_CODE_MISMATCH_K_SIZES = 101,
  SOURMASH_ERROR_CODE_MISMATCH_DNA_PROT = 102,
  SOURMASH_ERROR_CODE_MISMATCH_SCALED = 103,
  SOURMASH_ERROR_CODE_MISMATCH_SEED = 104,
  SOURMASH_ERROR_CODE_MISMATCH_SIGNATURE_TYPE = 105,
  SOURMASH_ERROR_CODE_NON_EMPTY_MIN_HASH = 106,
  SOURMASH_ERROR_CODE_MISMATCH_NUM = 107,
  SOURMASH_ERROR_CODE_INVALID_DNA = 1101,
  SOURMASH_ERROR_CODE_INVALID_PROT = 1102,
  SOURMASH_ERROR_CODE_INVALID_CODON_LENGTH = 1103,
  SOURMASH_ERROR_CODE_INVALID_HASH_FUNCTION = 1104,
  SOURMASH_ERROR_CODE_READ_DATA = 1201,
  SOURMASH_ERROR_CODE_STORAGE = 1202,
  SOURMASH_ERROR_CODE_HLL_PRECISION_BOUNDS = 1301,
  SOURMASH_ERROR_CODE_IO = 100001,
  SOURMASH_ERROR_CODE_UTF8_ERROR = 100002,
  SOURMASH_ERROR_CODE_PARSE_INT = 100003,
  SOURMASH_ERROR_CODE_SERDE_ERROR = 100004,
  SOURMASH_ERROR_CODE_NIFFLER_ERROR = 100005,
};
typedef uint32_t SourmashErrorCode;

typedef struct SourmashComputeParameters SourmashComputeParameters;

typedef struct SourmashHyperLogLog SourmashHyperLogLog;

typedef struct SourmashKmerMinHash SourmashKmerMinHash;

typedef struct SourmashNodegraph SourmashNodegraph;

typedef struct SourmashSignature SourmashSignature;

/**
 * Represents a string.
 */
typedef struct {
  /**
   * Pointer to the UTF-8 encoded string data.
   */
  char *data;
  /**
   * The length of the string pointed to by `data`.
   */
  uintptr_t len;
  /**
   * Indicates that the string is owned and must be freed.
   */
  bool owned;
} SourmashStr;

bool computeparams_dayhoff(const SourmashComputeParameters *ptr);

bool computeparams_dna(const SourmashComputeParameters *ptr);

void computeparams_free(SourmashComputeParameters *ptr);

bool computeparams_hp(const SourmashComputeParameters *ptr);

const uint32_t *computeparams_ksizes(const SourmashComputeParameters *ptr, uintptr_t *size);

void computeparams_ksizes_free(uint32_t *ptr, uintptr_t insize);

SourmashComputeParameters *computeparams_new(void);

uint32_t computeparams_num_hashes(const SourmashComputeParameters *ptr);

bool computeparams_protein(const SourmashComputeParameters *ptr);

uint64_t computeparams_scaled(const SourmashComputeParameters *ptr);

uint64_t computeparams_seed(const SourmashComputeParameters *ptr);

void computeparams_set_dayhoff(SourmashComputeParameters *ptr, bool v);

void computeparams_set_dna(SourmashComputeParameters *ptr, bool v);

void computeparams_set_hp(SourmashComputeParameters *ptr, bool v);

void computeparams_set_ksizes(SourmashComputeParameters *ptr,
                              const uint32_t *ksizes_ptr,
                              uintptr_t insize);

void computeparams_set_num_hashes(SourmashComputeParameters *ptr, uint32_t num);

void computeparams_set_protein(SourmashComputeParameters *ptr, bool v);

void computeparams_set_scaled(SourmashComputeParameters *ptr, uint64_t scaled);

void computeparams_set_seed(SourmashComputeParameters *ptr, uint64_t new_seed);

void computeparams_set_track_abundance(SourmashComputeParameters *ptr, bool v);

bool computeparams_track_abundance(const SourmashComputeParameters *ptr);

uint64_t hash_murmur(const char *kmer, uint64_t seed);

void hll_add_hash(SourmashHyperLogLog *ptr, uint64_t hash);

void hll_add_sequence(SourmashHyperLogLog *ptr, const char *sequence, uintptr_t insize, bool force);

uintptr_t hll_cardinality(const SourmashHyperLogLog *ptr);

double hll_containment(const SourmashHyperLogLog *ptr, const SourmashHyperLogLog *optr);

void hll_free(SourmashHyperLogLog *ptr);

SourmashHyperLogLog *hll_from_buffer(const char *ptr, uintptr_t insize);

SourmashHyperLogLog *hll_from_path(const char *filename);

uintptr_t hll_intersection_size(const SourmashHyperLogLog *ptr, const SourmashHyperLogLog *optr);

uintptr_t hll_ksize(const SourmashHyperLogLog *ptr);

uintptr_t hll_matches(const SourmashHyperLogLog *ptr, const SourmashKmerMinHash *mh_ptr);

void hll_merge(SourmashHyperLogLog *ptr, const SourmashHyperLogLog *optr);

SourmashHyperLogLog *hll_new(void);

void hll_save(const SourmashHyperLogLog *ptr, const char *filename);

double hll_similarity(const SourmashHyperLogLog *ptr, const SourmashHyperLogLog *optr);

const uint8_t *hll_to_buffer(const SourmashHyperLogLog *ptr, uintptr_t *size);

void hll_update_mh(SourmashHyperLogLog *ptr, const SourmashKmerMinHash *optr);

SourmashHyperLogLog *hll_with_error_rate(double error_rate, uintptr_t ksize);

void kmerminhash_add_from(SourmashKmerMinHash *ptr, const SourmashKmerMinHash *other);

void kmerminhash_add_hash(SourmashKmerMinHash *ptr, uint64_t h);

void kmerminhash_add_hash_with_abundance(SourmashKmerMinHash *ptr, uint64_t h, uint64_t abundance);

void kmerminhash_add_many(SourmashKmerMinHash *ptr, const uint64_t *hashes_ptr, uintptr_t insize);

void kmerminhash_add_protein(SourmashKmerMinHash *ptr, const char *sequence);

void kmerminhash_add_sequence(SourmashKmerMinHash *ptr, const char *sequence, bool force);

void kmerminhash_add_word(SourmashKmerMinHash *ptr, const char *word);

double kmerminhash_angular_similarity(const SourmashKmerMinHash *ptr,
                                      const SourmashKmerMinHash *other);

void kmerminhash_clear(SourmashKmerMinHash *ptr);

uint64_t kmerminhash_count_common(const SourmashKmerMinHash *ptr,
                                  const SourmashKmerMinHash *other,
                                  bool downsample);

bool kmerminhash_dayhoff(const SourmashKmerMinHash *ptr);

void kmerminhash_disable_abundance(SourmashKmerMinHash *ptr);

void kmerminhash_enable_abundance(SourmashKmerMinHash *ptr);

void kmerminhash_free(SourmashKmerMinHash *ptr);

const uint64_t *kmerminhash_get_abunds(SourmashKmerMinHash *ptr, uintptr_t *size);

const uint64_t *kmerminhash_get_mins(const SourmashKmerMinHash *ptr, uintptr_t *size);

uintptr_t kmerminhash_get_mins_size(const SourmashKmerMinHash *ptr);

HashFunctions kmerminhash_hash_function(const SourmashKmerMinHash *ptr);

void kmerminhash_hash_function_set(SourmashKmerMinHash *ptr, HashFunctions hash_function);

bool kmerminhash_hp(const SourmashKmerMinHash *ptr);

uint64_t kmerminhash_intersection(const SourmashKmerMinHash *ptr, const SourmashKmerMinHash *other);

bool kmerminhash_is_compatible(const SourmashKmerMinHash *ptr, const SourmashKmerMinHash *other);

bool kmerminhash_is_protein(const SourmashKmerMinHash *ptr);

double kmerminhash_jaccard(const SourmashKmerMinHash *ptr, const SourmashKmerMinHash *other);

uint32_t kmerminhash_ksize(const SourmashKmerMinHash *ptr);

uint64_t kmerminhash_max_hash(const SourmashKmerMinHash *ptr);

SourmashStr kmerminhash_md5sum(const SourmashKmerMinHash *ptr);

void kmerminhash_merge(SourmashKmerMinHash *ptr, const SourmashKmerMinHash *other);

SourmashKmerMinHash *kmerminhash_new(uint64_t scaled,
                                     uint32_t k,
                                     HashFunctions hash_function,
                                     uint64_t seed,
                                     bool track_abundance,
                                     uint32_t n);

uint32_t kmerminhash_num(const SourmashKmerMinHash *ptr);

void kmerminhash_remove_hash(SourmashKmerMinHash *ptr, uint64_t h);

void kmerminhash_remove_many(SourmashKmerMinHash *ptr,
                             const uint64_t *hashes_ptr,
                             uintptr_t insize);

uint64_t kmerminhash_seed(const SourmashKmerMinHash *ptr);

void kmerminhash_set_abundances(SourmashKmerMinHash *ptr,
                                const uint64_t *hashes_ptr,
                                const uint64_t *abunds_ptr,
                                uintptr_t insize,
                                bool clear);

double kmerminhash_similarity(const SourmashKmerMinHash *ptr,
                              const SourmashKmerMinHash *other,
                              bool ignore_abundance,
                              bool downsample);

void kmerminhash_slice_free(uint64_t *ptr, uintptr_t insize);

bool kmerminhash_track_abundance(const SourmashKmerMinHash *ptr);

double nodegraph_angular_similarity_upper_bound(const SourmashNodegraph *ptr,
                                                const SourmashKmerMinHash *mh_ptr);

void nodegraph_buffer_free(uint8_t *ptr, uintptr_t insize);

bool nodegraph_count(SourmashNodegraph *ptr, uint64_t h);

bool nodegraph_count_kmer(SourmashNodegraph *ptr, const char *kmer);

double nodegraph_expected_collisions(const SourmashNodegraph *ptr);

void nodegraph_free(SourmashNodegraph *ptr);

SourmashNodegraph *nodegraph_from_buffer(const char *ptr, uintptr_t insize);

SourmashNodegraph *nodegraph_from_path(const char *filename);

uintptr_t nodegraph_get(const SourmashNodegraph *ptr, uint64_t h);

uintptr_t nodegraph_get_kmer(const SourmashNodegraph *ptr, const char *kmer);

const uint64_t *nodegraph_hashsizes(const SourmashNodegraph *ptr, uintptr_t *size);

uintptr_t nodegraph_ksize(const SourmashNodegraph *ptr);

uintptr_t nodegraph_matches(const SourmashNodegraph *ptr, const SourmashKmerMinHash *mh_ptr);

SourmashNodegraph *nodegraph_new(void);

uintptr_t nodegraph_noccupied(const SourmashNodegraph *ptr);

uintptr_t nodegraph_ntables(const SourmashNodegraph *ptr);

void nodegraph_save(const SourmashNodegraph *ptr, const char *filename);

const uint8_t *nodegraph_to_buffer(const SourmashNodegraph *ptr,
                                   uint8_t compression,
                                   uintptr_t *size);

void nodegraph_update(SourmashNodegraph *ptr, const SourmashNodegraph *optr);

void nodegraph_update_mh(SourmashNodegraph *ptr, const SourmashKmerMinHash *optr);

SourmashNodegraph *nodegraph_with_tables(uintptr_t ksize,
                                         uintptr_t starting_size,
                                         uintptr_t n_tables);

void signature_add_protein(SourmashSignature *ptr, const char *sequence);

void signature_add_sequence(SourmashSignature *ptr, const char *sequence, bool force);

bool signature_eq(const SourmashSignature *ptr, const SourmashSignature *other);

SourmashKmerMinHash *signature_first_mh(const SourmashSignature *ptr);

void signature_free(SourmashSignature *ptr);

SourmashSignature *signature_from_params(const SourmashComputeParameters *ptr);

SourmashStr signature_get_filename(const SourmashSignature *ptr);

SourmashStr signature_get_license(const SourmashSignature *ptr);

SourmashKmerMinHash **signature_get_mhs(const SourmashSignature *ptr, uintptr_t *size);

SourmashStr signature_get_name(const SourmashSignature *ptr);

uintptr_t signature_len(const SourmashSignature *ptr);

SourmashSignature *signature_new(void);

void signature_push_mh(SourmashSignature *ptr, const SourmashKmerMinHash *other);

SourmashStr signature_save_json(const SourmashSignature *ptr);

void signature_set_filename(SourmashSignature *ptr, const char *name);

void signature_set_mh(SourmashSignature *ptr, const SourmashKmerMinHash *other);

void signature_set_name(SourmashSignature *ptr, const char *name);

SourmashSignature **signatures_load_buffer(const char *ptr,
                                           uintptr_t insize,
                                           bool _ignore_md5sum,
                                           uintptr_t ksize,
                                           const char *select_moltype,
                                           uintptr_t *size);

SourmashSignature **signatures_load_path(const char *ptr,
                                         bool _ignore_md5sum,
                                         uintptr_t ksize,
                                         const char *select_moltype,
                                         uintptr_t *size);

const uint8_t *signatures_save_buffer(const SourmashSignature *const *ptr,
                                      uintptr_t size,
                                      uint8_t compression,
                                      uintptr_t *osize);

char sourmash_aa_to_dayhoff(char aa);

char sourmash_aa_to_hp(char aa);

/**
 * Clears the last error.
 */
void sourmash_err_clear(void);

/**
 * Returns the panic information as string.
 */
SourmashStr sourmash_err_get_backtrace(void);

/**
 * Returns the last error code.
 *
 * If there is no error, 0 is returned.
 */
SourmashErrorCode sourmash_err_get_last_code(void);

/**
 * Returns the last error message.
 *
 * If there is no error an empty string is returned.  This allocates new memory
 * that needs to be freed with `sourmash_str_free`.
 */
SourmashStr sourmash_err_get_last_message(void);

/**
 * Initializes the library
 */
void sourmash_init(void);

/**
 * Frees a sourmash str.
 *
 * If the string is marked as not owned then this function does not
 * do anything.
 */
void sourmash_str_free(SourmashStr *s);

/**
 * Creates a sourmash str from a c string.
 *
 * This sets the string to owned.  In case it's not owned you either have
 * to make sure you are not freeing the memory or you need to set the
 * owned flag to false.
 */
SourmashStr sourmash_str_from_cstr(const char *s);

char sourmash_translate_codon(const char *codon);

#endif /* SOURMASH_H_INCLUDED */
