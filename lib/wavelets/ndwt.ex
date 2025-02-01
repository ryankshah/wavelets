defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform
  """

  alias Wavelets.{DWT, Filter, IDWT}

  @doc """
  Performs n-dimensional forward discrete wavelet transform
  """
  @spec forward(
          list() | list(list()) | list(list(list())),
          Filter.t(),
          pos_integer(),
          Wavelets.precision()
        ) ::
          {list(), list()}
  def forward(signal, filter, dims, _precision \\ :double) do
    case dims do
      1 when is_list(signal) ->
        DWT.forward_1d(signal, filter)

      2 when is_list(hd(signal)) ->
        DWT.forward_2d(signal, filter)

      n when n > 2 ->
        do_forward_nd(signal, filter, dims)

      _ ->
        raise ArgumentError, "Invalid signal dimensions or structure"
    end
  end

  @doc """
  Performs n-dimensional inverse discrete wavelet transform
  """
  @spec inverse(list(), list(), Filter.t(), pos_integer(), Wavelets.precision()) ::
          list()
  def inverse(approximation, details, filter, dims, _precision \\ :double) do
    case dims do
      1 ->
        IDWT.inverse_1d(approximation, details, filter)

      2 ->
        IDWT.inverse_2d(approximation, details, filter)

      n when n > 2 ->
        do_inverse_nd(approximation, details, filter, dims)

      _ ->
        raise ArgumentError, "Invalid dimensions"
    end
  end

  # Helper functions
  defp do_forward_nd(signal, filter, dims) when dims > 2 do
    # Apply transform along first dimension
    {approx_temp, details_temp} =
      Enum.map(signal, fn slice ->
        if is_list(slice) do
          forward(slice, filter, dims - 1)
        else
          raise ArgumentError,
                "Invalid signal structure for n-dimensional transform"
        end
      end)
      |> Enum.unzip()

    {approx_temp, details_temp}
  end

  defp do_inverse_nd(approximation, details, filter, dims) when dims > 2 do
    Enum.zip(approximation, details)
    |> Enum.map(fn {approx_slice, detail_slice} ->
      inverse(approx_slice, detail_slice, filter, dims - 1)
    end)
  end
end
